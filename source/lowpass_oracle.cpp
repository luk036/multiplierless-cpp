#include <cmath>
#include <multiplierless/lowpass_oracle.hpp>
#include <numbers>
#include <optional>

#ifndef M_PI
constexpr double M_PI = std::numbers::pi;
#endif

namespace {

    /// @brief Dot product of one constraint row with the variable vector x.
    /// @param[in] mat Constraint matrix
    /// @param[in] row Row index
    /// @param[in] x   Variable vector
    /// @return Row dot product
    auto dot_row(const Arr& mat, std::size_t row, const Arr& x) -> double {
        const double* row_data = mat.data() + row * mat.cols();
        const double* x_data = x.data();
        const auto n = x.size();
        double sum = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            sum += row_data[j] * x_data[j];
        }
        return sum;
    }

    /// @brief Template-Method skeleton: scan the rows of `mat` in round-robin
    /// order and return the first violating cut reported by `check`, or nullopt.
    ///
    /// The `check` callback receives (row, dot) and returns an optional cut;
    /// returning nullopt continues the scan.
    template <typename Check> auto scan_constraints(const Arr& mat, RoundRobin& rr, const Arr& x,
                                                    Check&& check) -> std::optional<ParallelCut> {
        const auto n = mat.rows();
        for (auto i = 0U; i != n; ++i) {
            const auto k = rr.next();
            if (auto cut = check(k, dot_row(mat, k, x))) {
                return cut;
            }
        }
        return std::nullopt;
    }

}  // namespace

/**
 * @brief Default constructor using built-in default filter specs.
 *
 * Delegates to the full constructor with defaults:
 *   - Normalized passband edge:  0.12
 *   - Normalized stopband edge:  0.20
 *   - Passband ripple:           0.125
 *   - Stopband attenuation:      0.125
 *   - Discretization factor:     15
 *
 * @param[in] argN Filter order (number of FIR coefficients).
 */
filter_design_construct::filter_design_construct(int argN)
    : filter_design_construct(argN, 0.12, 0.20, 0.125, 0.125, 15) {}

/**
 * @brief Construct a filter design with explicit specifications.
 *
 * Builds the constraint matrices Ap (passband), As (stopband), Anr
 * (non-redundant) and the squared bounds Lpsq, Upsq, Spsq from the
 * given analog filter parameters.
 *
 * @param[in] argN             Filter order.
 * @param[in] wpass_norm       Normalized passband edge frequency (×π rad/sample).
 * @param[in] wstop_norm       Normalized stopband edge frequency (×π rad/sample).
 * @param[in] passband_ripple  Allowed passband ripple (linear).
 * @param[in] stopband_attn    Stopband attenuation (linear).
 * @param[in] discretization_factor  Grid density multiplier (m = factor × N).
 */
filter_design_construct::filter_design_construct(int argN, double wpass_norm, double wstop_norm,
                                                 double passband_ripple, double stopband_attn,
                                                 int discretization_factor)
    : N(argN) {
    const auto wpass = wpass_norm * M_PI;
    const auto wstop = wstop_norm * M_PI;
    const auto delta = 20 * std::log10(1 + passband_ripple);
    const auto delta2 = 20 * std::log10(stopband_attn);
    const auto m = discretization_factor * this->N;
    const auto m_sz = static_cast<size_t>(m);
    const auto N1 = static_cast<size_t>(this->N - 1);
    const auto w = linspace(0, M_PI, m_sz);  // omega
    // passband 0 <= w <= w_pass
    const auto Lp = std::pow(10, -delta / 20);
    const auto Up = std::pow(10, +delta / 20);
    // A is the matrix used to compute the power spectrum
    // A(w,:) = [1 2*cos(w) 2*cos(2*w) ... 2*cos((this->N-1)*w)]
    Arr An = zeros(m_sz, N1);
    for (size_t i = 0; i != m_sz; ++i) {
        for (size_t j = 0; j != N1; ++j) {
            An(i, j) = 2.0 * std::cos(w(i) * static_cast<double>(j + 1));
        }
    }
    Arr A = concatenate(ones(m_sz, 1), An, 1);
    const auto ind_p = where(w <= wpass)[0];  // passband
    auto ind_p_size = ind_p.size();
    this->Ap = view(A, Range(0, ind_p_size), Range(Range::ALL));
    // stopband (w >= w_stop)
    const auto ind_s = where(w >= wstop)[0];  // stopband
    const auto Sp = std::pow(10, delta2 / 20);
    auto ind_s_0 = static_cast<size_t>(ind_s[0]);
    this->As = view(A, Range(ind_s_0, Range::ALL), Range(Range::ALL));
    // Remove redundant contraints
    auto ind_p_last = static_cast<size_t>(ind_p[ind_p_size - 1]);
    auto ind_end = static_cast<size_t>(ind_s[0]);
    this->Anr = view(A, Range(ind_p_last + 1, ind_end), Range(Range::ALL));
    this->Lpsq = Lp * Lp;
    this->Upsq = Up * Up;
    this->Spsq = Sp * Sp;
}

/**
 * @brief Assess the optimization for the given autocorrelation coefficients.
 *
 * Evaluates the non-negative-real constraint, passband constraints (Upsq/Lpsq),
 * and stopband constraint (Spsq) using a round-robin traversal of the constraint
 * matrices. Returns the cutting-plane (gradient + objective) when a constraint is
 * violated, or signals optimality when all constraints are satisfied.
 *
 * @param[in] x A 1-dimensional array representing the optimization variables.
 * @param[in,out] Spsq On input, the target stopband attenuation squared.
 *                      On output, the achieved maximum stopband value.
 *
 * @return The function `assess_optim` returns a tuple containing a `ParallelCut` object and a
 * boolean value.
 */
auto LowpassOracle::assess_optim(const Arr& x, double& Spsq) -> std::tuple<ParallelCut, bool> {
    // 1.0 nonnegative-real constraint
    // case 1,
    if (x(0) < 0) {
        Arr g = zeros(x.size());
        g(0) = -1.;
        auto f = Vec{-x(0)};
        return {{std::move(g), std::move(f)}, false};
    }

    // 2.0 passband constraints
    if (auto cut = scan_constraints(this->_Fdc.Ap, this->_rr_ap, x,
                                    [&](std::size_t k, double v) -> std::optional<ParallelCut> {
                                        if (v > this->_Fdc.Upsq) {
                                            // Calculate: f = v - Upsq;
                                            Arr g(this->_Fdc.Ap.cols());
                                            for (std::size_t j = 0; j < this->_Fdc.Ap.cols(); ++j) {
                                                g(j) = this->_Fdc.Ap(k, j);
                                            }
                                            Vec f{v - this->_Fdc.Upsq, v - this->_Fdc.Lpsq};
                                            return ParallelCut{std::move(g), std::move(f)};
                                        }
                                        if (v < this->_Fdc.Lpsq) {
                                            // Calculate: f = Lpsq - v;
                                            Arr g(this->_Fdc.Ap.cols());
                                            for (std::size_t j = 0; j < this->_Fdc.Ap.cols(); ++j) {
                                                g(j) = -this->_Fdc.Ap(k, j);
                                            }
                                            Vec f{-v + this->_Fdc.Lpsq, -v + this->_Fdc.Upsq};
                                            return ParallelCut{std::move(g), std::move(f)};
                                        }
                                        return std::nullopt;
                                    })) {
        return {std::move(*cut), false};
    }

    // 3.0 stopband constraint
    auto fmax = -1.e100;
    std::size_t imax = 0U;
    if (auto cut = scan_constraints(this->_Fdc.As, this->_rr_as, x,
                                    [&](std::size_t k, double v) -> std::optional<ParallelCut> {
                                        if (v > Spsq) {
                                            // Calculate: f = v - Spsq
                                            Arr g(this->_Fdc.As.cols());
                                            for (std::size_t j = 0; j < this->_Fdc.As.cols(); ++j) {
                                                g(j) = this->_Fdc.As(k, j);
                                            }
                                            Vec f{v - Spsq, v};
                                            return ParallelCut{std::move(g), std::move(f)};
                                        }
                                        if (v < 0) {
                                            // Calculate: f = v - Spsq
                                            Arr g(this->_Fdc.As.cols());
                                            for (std::size_t j = 0; j < this->_Fdc.As.cols(); ++j) {
                                                g(j) = -this->_Fdc.As(k, j);
                                            }
                                            Vec f{-v, -v + Spsq};
                                            return ParallelCut{std::move(g), std::move(f)};
                                        }
                                        if (v > fmax) {
                                            fmax = v;
                                            imax = k;
                                        }
                                        return std::nullopt;
                                    })) {
        return {std::move(*cut), false};
    }

    // 4.0 nonnegative-real constraint on the non-redundant rows
    if (auto cut = scan_constraints(this->_Fdc.Anr, this->_rr_anr, x,
                                    [&](std::size_t k, double v) -> std::optional<ParallelCut> {
                                        if (v < 0.0) {
                                            Vec f{-v};
                                            Arr g(this->_Fdc.Anr.cols());
                                            for (std::size_t j = 0; j < this->_Fdc.Anr.cols();
                                                 ++j) {
                                                g(j) = -this->_Fdc.Anr(k, j);
                                            }
                                            return ParallelCut{std::move(g), std::move(f)};
                                        }
                                        return std::nullopt;
                                    })) {
        return {std::move(*cut), false};
    }

    // Begin objective function
    Spsq = fmax;
    Vec f{0.0, fmax};
    Arr g(this->_Fdc.As.cols());
    for (std::size_t j = 0; j < this->_Fdc.As.cols(); ++j) {
        g(j) = this->_Fdc.As(imax, j);
    }
    return {{std::move(g), std::move(f)}, true};
}
