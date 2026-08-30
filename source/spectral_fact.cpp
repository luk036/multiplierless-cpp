#include <algorithm>
#include <cmath>
#include <complex>
#include <ellalgo/arr.hpp>
#include <ginger/aberth.hpp>
#include <ginger/config.hpp>
#include <multiplierless/fftw_helper.hpp>
#include <multiplierless/spectral_fact.hpp>
#include <numbers>
#include <utility>
#include <vector>

#ifndef M_PI
constexpr double M_PI = std::numbers::pi;
#endif

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
    std::ranges::reverse(coeffs);

    auto zs = initial_aberth_autocorr(coeffs);
    ginger::Options opts;
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
 * @brief Spectral factorization via FFT / Hilbert transform.
 *
 * Over-samples the frequency response by a factor of 100, computes
 * \f$ \alpha = \frac{1}{2}\ln|R(\omega)| \f$, applies the Hilbert
 * transform to obtain the minimum-phase log-magnitude / phase pair,
 * and returns the inverse FFT of \f$ e^{\alpha + j\phi} \f$.
 *
 * Uses FFT (via rfft) to compute the frequency response instead of a stored matrix.
 *
 * @param[in] r Autocorrelation sequence (top half).
 * @return Minimum-phase impulse response h.
 */
auto spectral_fact_fft(const Arr& r) -> Arr {
    const auto n = static_cast<int>(r.size());
    const auto mult_factor = 100;
    const auto m = mult_factor * n;
    const auto m_sz = static_cast<size_t>(m);

    // Compute R(ω) = r₀ + 2·Σ_{k=1}^{n-1} r[k]·cos(k·ω) via FFT instead of matrix multiply.
    // Zero-pad r to length m and take rfft: S[i] = Σ r[k]·exp(-j·k·ω_i)
    // Then R[i] = 2·Re(S[i]) - r₀, since Re(S[i]) = r₀ + Σ r[k]·cos(k·ω_i)
    auto pad = zeros(m_sz);
    for (size_t i = 0; i < r.size(); ++i) pad(i) = r(i);
    auto S = rfft(pad);

    Arr R(m_sz);
    const double r0 = r(0);
    const auto half = m_sz / 2;
    for (size_t i = 0; i <= half; ++i) {
        R(i) = 2.0 * S[i].real() - r0;
    }
    for (size_t i = half + 1; i < m_sz; ++i) {
        R(i) = R(m_sz - i);
    }

    auto min_val = *std::ranges::min_element(R);
    if (min_val <= 0) {
        for (size_t i = 0; i < R.size(); ++i)
            if (R(i) <= 0) R(i) = 1e-10;
    }

    Arr alpha = 0.5 * log(abs(R));
    auto alphatmp = fft(cast_to_complex(alpha));
    auto ind = static_cast<size_t>(m) / 2;
    for (auto i = ind; std::cmp_less(i, m); ++i) alphatmp[i] = -alphatmp[i];
    alphatmp[0] = {0.0, 0.0};
    alphatmp[ind] = {0.0, 0.0};

    const std::complex<double> j_{0, 1};
    Arr phi = ifft(j_ * alphatmp);
    auto alpha1 = view(alpha, Range(0, m_sz, static_cast<size_t>(mult_factor)));
    auto phi1 = view(phi, Range(0, m_sz, static_cast<size_t>(mult_factor)));
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
