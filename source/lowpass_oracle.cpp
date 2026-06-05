#include <cmath>
#include <multiplierless/lowpass_oracle.hpp>

using Vec = std::valarray<double>;
using ParallelCut = std::pair<Arr, Vec>;

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846264338327950288;
#endif

/**
 * The above function is a constructor for a filter design class that initializes various parameters
 * and matrices used in the filter design process.
 *
 * @param[in] argN The parameter `argN` represents the value of N, which is the order of the filter.
 * It determines the number of filter coefficients and the complexity of the filter design.
 */
filter_design_construct::filter_design_construct(int argN)
    : filter_design_construct(argN, 0.12, 0.20, 0.125, 0.125, 15) {}

filter_design_construct::filter_design_construct(int argN, double wpass_norm, double wstop_norm,
                                                  double passband_ripple, double stopband_attn,
                                                  int discretization_factor)
    : N(argN) {
    const auto wpass = wpass_norm * M_PI;
    const auto wstop = wstop_norm * M_PI;
    const auto delta = 20 * std::log10(1 + passband_ripple);
    const auto delta2 = 20 * std::log10(stopband_attn);
    const auto m = discretization_factor * this->N;
    const auto w = linspace(0, M_PI, size_t(m));  // omega
    // passband 0 <= w <= w_pass
    const auto Lp = std::pow(10, -delta / 20);
    const auto Up = std::pow(10, +delta / 20);
    // A is the matrix used to compute the power spectrum
    // A(w,:) = [1 2*cos(w) 2*cos(2*w) ... 2*cos((this->N-1)*w)]
    Arr An = zeros(m, this->N - 1);
    for (auto i = 0; i != m; ++i) {
        for (auto j = 0; j != this->N - 1; ++j) {
            An(i, j) = 2.0 * std::cos(w(i) * (j + 1));
        }
    }
    Arr A = concatenate(ones(m, 1), An, 1);
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
 * The function assess_optim in the LowpassOracle class assesses the optimization of a given input
 * vector x based on various constraints and returns a tuple containing the gradient and objective
 * function values, along with a boolean indicating whether the optimization is complete.
 *
 * @param[in] x A 1-dimensional array representing the optimization variables.
 * @param[in] Spsq Spsq is a reference to a double variable. It is used to store the maximum value
 * of the stopband constraint.
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

    // case 2,
    // 2.0 passband constraints
    auto N = this->_Fdc.Ap.rows();

    auto dot_row = [&](const Arr& mat, size_t row) -> double {
        double sum = 0.0;
        for (size_t j = 0; j < x.size(); ++j) {
            sum += mat(row, j) * x(j);
        }
        return sum;
    };

    auto k = this->_i_Ap;
    for (auto i = 0U; i != N; ++i, ++k) {
        if (k == N) {
            k = 0;  // round robin
        }
        auto v = dot_row(this->_Fdc.Ap, k);
        if (v > this->_Fdc.Upsq) {
            // Calculate: f = v - Upsq;
            Arr g(this->_Fdc.Ap.cols());
            for (size_t j = 0; j < this->_Fdc.Ap.cols(); ++j) g(j) = this->_Fdc.Ap(k, j);
            Vec f{v - this->_Fdc.Upsq, v - this->_Fdc.Lpsq};
            this->_i_Ap = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
        if (v < this->_Fdc.Lpsq) {
            // Calculate: f = Lpsq - v;
            Arr g(this->_Fdc.Ap.cols());
            for (size_t j = 0; j < this->_Fdc.Ap.cols(); ++j) g(j) = -this->_Fdc.Ap(k, j);
            Vec f{-v + this->_Fdc.Lpsq, -v + this->_Fdc.Upsq};
            this->_i_Ap = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
    }

    // case 3,
    // 3.0 stopband constraint
    N = this->_Fdc.As.rows();
    auto fmax = -1.e100;
    size_t imax = 0U;
    k = this->_i_As;
    for (auto i = 0U; i != N; ++i, ++k) {
        if (k == N) {
            k = 0;  // round robin
        }
        auto v = dot_row(this->_Fdc.As, k);
        if (v > Spsq) {
            // Calculate: f = v - Spsq
            Arr g(this->_Fdc.As.cols());
            for (size_t j = 0; j < this->_Fdc.As.cols(); ++j) g(j) = this->_Fdc.As(k, j);
            Vec f{v - Spsq, v};
            this->_i_As = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
        if (v < 0) {
            // Calculate: f = v - Spsq
            Arr g(this->_Fdc.As.cols());
            for (size_t j = 0; j < this->_Fdc.As.cols(); ++j) g(j) = -this->_Fdc.As(k, j);
            Vec f{-v, -v + Spsq};
            this->_i_As = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
        if (v > fmax) {
            fmax = v;
            imax = k;
        }
    }

    // case 4,
    // 1.0 nonnegative-real constraint
    N = this->_Fdc.Anr.rows();
    k = this->_i_Anr;
    for (auto i = 0U; i != N; ++i, ++k) {
        if (k == N) {
            k = 0;  // round robin
        }
        auto v = dot_row(this->_Fdc.Anr, k);
        if (v < 0.0) {
            Vec f{-v};
            Arr g(this->_Fdc.Anr.cols());
            for (size_t j = 0; j < this->_Fdc.Anr.cols(); ++j) g(j) = -this->_Fdc.Anr(k, j);
            this->_i_Anr = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
    }

    // Begin objective function
    Spsq = fmax;
    Vec f{0.0, fmax};
    Arr g(this->_Fdc.As.cols());
    for (size_t j = 0; j < this->_Fdc.As.cols(); ++j) g(j) = this->_Fdc.As(imax, j);
    return {{std::move(g), std::move(f)}, true};
}
