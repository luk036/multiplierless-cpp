#include <algorithm>
#include <cmath>
#include <complex>
#include <ellalgo/arr.hpp>
#include <ginger/aberth.hpp>
#include <ginger/config.hpp>
#include <multiplierless/fftw_helper.hpp>
#include <vector>

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846264338327950288;
#endif

auto spectral_fact_fft(const Arr& r) -> Arr;

/**
 * @brief Spectral factorization via root-finding.
 *
 * Builds a symmetric polynomial from the autocorrelation coefficients,
 * finds its roots using the Aberth-Ehrlich method (via ginger library),
 * selects roots inside the unit circle, and reconstructs the minimum-phase
 * impulse response. Normalises the energy to match r(0).
 *
 * @param[in] r         Autocorrelation sequence (top half).
 * @param[in] tolerance Convergence tolerance for the Aberth solver.
 * @return The minimum-phase impulse response coefficients h.
 */
auto spectral_fact_root(const Arr& r, double tolerance) -> Arr {
    const auto n = r.size();
    const auto deg = 2 * n - 2;
    std::vector<double> coeffs(deg + 1, 0.0);
    coeffs[0] = r(n - 1);
    for (size_t i = 0; i < n - 1; ++i) coeffs[i + 1] = 2.0 * r(n - 2 - i);
    for (size_t i = 0; i < n - 2; ++i) coeffs[deg - i - 1] = 2.0 * r(n - 2 - i);
    coeffs[n - 1] = 2.0 * r(0);
    coeffs[deg] = r(n - 1);
    std::reverse(coeffs.begin(), coeffs.end());

    auto zs = initial_aberth_autocorr(coeffs);
    Options opts;
    opts.tolerance = tolerance;
    opts.max_iters = 500;
    aberth_autocorr(coeffs, zs, opts);

    std::vector<std::complex<double>> inside;
    for (auto& z : zs) {
        if (std::abs(z) < 1.0)
            inside.push_back(z);
        else
            inside.push_back(1.0 / z);
    }

    auto hc = poly_from_roots(inside);
    double eh = 0.0;
    for (auto c : hc) eh += c * c;
    const auto norm = std::sqrt(r(0) / eh);
    for (auto& c : hc) c *= norm;

    Arr h(n);
    for (size_t i = 0; i < n && i < hc.size(); ++i) h(i) = hc[i];
    return h;
}

/**
 * @brief Spectral factorization (convenience wrapper).
 *
 * Delegates to spectral_fact_fft.
 * @param[in] r Autocorrelation sequence.
 * @return Minimum-phase impulse response h.
 */
auto spectral_fact(const Arr& r) -> Arr { return spectral_fact_fft(r); }

/**
 * @brief Spectral factorization via FFT / Hilbert transform.
 *
 * Over-samples the frequency response by a factor of 100, computes
 * \f$ \alpha = \frac{1}{2}\ln|R(\omega)| \f$, applies the Hilbert
 * transform to obtain the minimum-phase log-magnitude / phase pair,
 * and returns the inverse FFT of \f$ e^{\alpha + j\phi} \f$.
 *
 * Results are cached for repeated calls with the same filter order.
 *
 * @param[in] r Autocorrelation sequence (top half).
 * @return Minimum-phase impulse response h.
 */
auto spectral_fact_fft(const Arr& r) -> Arr {
    const auto n = static_cast<int>(r.size());
    const auto mult_factor = 100;
    const auto m = mult_factor * n;

    static int cached_n = 0;
    static Arr cached_A;
    if (n != cached_n) {
        const auto step_w = 2.0 * M_PI / static_cast<double>(m);
        Arr w(m);
        for (size_t i = 0; i < static_cast<size_t>(m); ++i) w(i) = static_cast<double>(i) * step_w;
        auto cols = arange(1.0, static_cast<double>(n));
        Arr An = 2.0 * cos(outer(w, cols));
        cached_A = concatenate(ones(m, 1), An, 1);
        cached_n = n;
    }
    const auto& A = cached_A;
    Arr R = dot(A, r);

    auto min_val = *std::min_element(R.begin(), R.end());
    if (min_val <= 0) {
        for (size_t i = 0; i < R.size(); ++i)
            if (R(i) <= 0) R(i) = 1e-10;
    }

    Arr alpha = 0.5 * log(abs(R));
    auto alphatmp = fft(cast_to_complex(alpha));
    auto ind = static_cast<size_t>(m) / 2;
    for (auto i = ind; i < m; ++i) alphatmp[i] = -alphatmp[i];
    alphatmp[0] = {0.0, 0.0};
    alphatmp[ind] = {0.0, 0.0};

    const std::complex<double> j_{0, 1};
    Arr phi = ifft(j_ * alphatmp);
    auto alpha1 = view(alpha, Range(0, m, static_cast<size_t>(mult_factor)));
    auto phi1 = view(phi, Range(0, m, static_cast<size_t>(mult_factor)));
    return ifft(exp(cast_to_complex(alpha1) + j_ * cast_to_complex(phi1)));
}

/**
 * @brief Inverse spectral factorization.
 *
 * Computes the autocorrelation sequence from an impulse response:
 * \f$ r(t) = \sum_{i=0}^{n-1-t} h(i+t) \cdot h(i) \f$.
 *
 * @param[in] h Impulse response coefficients.
 * @return Autocorrelation sequence r (same length as h).
 */
auto inverse_spectral_fact(const Arr& h) -> Arr {
    auto n = h.size();
    Arr r(n);
    for (size_t t = 0; t < n; ++t) {
        double sum = 0.0;
        for (size_t i = 0; i < n - t; ++i) sum += h(i + t) * h(i);
        r(t) = sum;
    }
    return r;
}
