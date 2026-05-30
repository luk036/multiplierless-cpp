#include <algorithm>
#include <cmath>
#include <complex>

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846264338327950288;
#endif

#include <ellalgo/arr.hpp>
#include <multiplierless/fftw_helper.hpp>

/* The `spectral_fact` function is performing spectral factorization using the
Kolmogorov 1939 approach. It takes an input vector `r` which represents the
top-half of the auto-correlation coefficients. It computes the minimum-phase
impulse response `h` that satisfies the given auto-correlation. The function
returns the impulse response `h` as an `Arr` object. */

/**
 * @brief Spectral factorization
 *
 * The spectral_fact function performs spectral factorization using the Kolmogorov 1939 approach to
 * compute the minimum-phase impulse response that satisfies a given auto-correlation.
 * (code follows pp. 232-233, Signal Analysis, by A. Papoulis)
 *
 * @param[in] r The parameter `r` is the top-half of the auto-correlation coefficients. It
 * represents the desired auto-correlation of the impulse response. It should be passed in as a
 * column vector.
 *
 * @return The function `spectral_fact` returns the impulse response `h` that gives the desired
 * auto-correlation.
 */
auto spectral_fact(const Arr& r) -> Arr {
    const auto n = int(r.size());
    const auto mult_factor = 100;
    const auto m = mult_factor * n;

    // Cache the cosine matrix A across calls (depends only on n, not on r)
    // Note: linspace with m points from 0 to 2π, endpoint excluded (Python convention)
    static int cached_n = 0;
    static Arr cached_A;
    if (n != cached_n) {
        double step_w = 2.0 * M_PI / static_cast<double>(m);
        Arr w(m);
        for (size_t i = 0; i < size_t(m); ++i) w(i) = static_cast<double>(i) * step_w;
        auto cols = arange(1.0, double(n));
        Arr An = 2.0 * cos(outer(w, cols));
        cached_A = concatenate(ones(m, 1), An, 1);
        cached_n = n;
    }
    const auto& A = cached_A;
    Arr R = dot(A, r);

    // Clamp near-zero freq response to avoid log(≤0)
    auto min_val = *std::min_element(R.begin(), R.end());
    if (min_val <= 0) {
        for (size_t i = 0; i < R.size(); ++i)
            if (R(i) <= 0) R(i) = 1e-10;
    }

    Arr alpha = 0.5 * log(abs(R));

    // Hilbert transform via full FFT (matches Python np.fft.fft / np.fft.ifft)
    auto alphatmp = fft(cast_to_complex(alpha));
    auto ind = size_t(m) / 2;
    // Negate the negative-frequency half (ind .. m-1)
    for (auto i = ind; i < m; ++i) alphatmp[i] = -alphatmp[i];
    alphatmp[0] = std::complex<double>(0.0, 0.0);
    alphatmp[ind] = std::complex<double>(0.0, 0.0);

    // phi = real(ifft(j * alphatmp))    — matches Python
    const std::complex<double> j_{0, 1};
    Arr phi = ifft(j_ * alphatmp);

    // Subsample alpha and phi by mult_factor
    auto alpha1 = view(alpha, Range(0, m, size_t(mult_factor)));
    auto phi1 = view(phi, Range(0, m, size_t(mult_factor)));

    // h = real(ifft(exp(alpha1 + j*phi1)))   — matches Python
    auto h_tmp = exp(cast_to_complex(alpha1) + j_ * cast_to_complex(phi1));
    Arr h = ifft(h_tmp);
    return h;
}

/**
 * The `inverse_spectral_fact` function takes an impulse response `h` as input
 * and computes the auto-correlation coefficients `r` that correspond to the
 * given impulse response. It returns the auto-correlation coefficients `r` as
 * an `Arr` object.
 *
 * Uses explicit O(n^2) convolution matching Python's np.convolve(h, h[::-1])[n-1:].
 *
 * @param[in] h The parameter `h` is the impulse response, which is a one-dimensional array or
 * vector representing the response of a system to an impulse input.
 *
 * @return The function `inverse_spectral_fact` returns the auto-correlation coefficients `r` as an
 * `Arr` object.
 */
auto inverse_spectral_fact(const Arr& h) -> Arr {
    auto n = h.size();
    Arr r(n);
    // r[t] = Σ_{i=0}^{n-1-t} h[i+t] * h[i]   for t = 0..n-1
    for (size_t t = 0; t < n; ++t) {
        double sum = 0.0;
        for (size_t i = 0; i < n - t; ++i) {
            sum += h(i + t) * h(i);
        }
        r(t) = sum;
    }
    return r;
}
