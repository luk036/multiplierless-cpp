// Disable svector on macOS to avoid Clang template ambiguity issues
// where long and unsigned long are both 64-bit
#ifdef __APPLE__
#    define XTENSOR_DISABLE_SVECTOR 1
#endif

#include <cmath>    // for cos, M_PI
#include <complex>  // for complex, operator*, operator-

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846264338327950288;
#endif

#include <xtensor-blas/xlinalg.hpp>  // for dot
#include <xtensor-fftw/basic.hpp>    // for irfft, rfft
#include <xtensor/xaccessible.hpp>   // for xaccessible
#include <xtensor/xarray.hpp>        // for xarray_container
#include <xtensor/xbroadcast.hpp>    // for xbroadcast
#include <xtensor/xbuilder.hpp>      // for zeros, concatenate, linspace
#include <xtensor/xcontainer.hpp>    // for xcontainer, xcontainer<>::...
#include <xtensor/xeval.hpp>         // for eval
#include <xtensor/xexception.hpp>    // for throw_concatenate_error
#include <xtensor/xfunction.hpp>     // for xfunction
#include <xtensor/xgenerator.hpp>    // for xgenerator
#include <xtensor/xiterator.hpp>     // for linear_begin
#include <xtensor/xlayout.hpp>       // for layout_type, layout_type::...
#include <xtensor/xmath.hpp>         // for abs, exp, log, sum, abs_fun
#include <xtensor/xoperation.hpp>    // for xfunction_type_t, operator*
#include <xtensor/xreducer.hpp>      // for xreducer
#include <xtensor/xslice.hpp>        // for range, xtuph, _
#include <xtensor/xutils.hpp>        // for accumulate
#include <xtensor/xview.hpp>         // for xview, view

using Arr = xt::xarray<double>;

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
    // length of the impulse response sequence
    const auto n = int(r.shape()[0]);

    // over-sampling factor
    const auto mult_factor = 20;  // should have mult_factor*(n) >> n
    const auto m = mult_factor * n;

    // Cache the cosine matrix A across calls (depends only on n, not on r)
    // This avoids O(m*n) recomputation on every optimization iteration.
    static int cached_n = 0;
    static Arr cached_A;
    if (n != cached_n) {
        // Use xt::linalg::outer for vectorized matrix construction
        Arr w = xt::linspace<double>(0, 2 * M_PI, size_t(m));
        auto cols = xt::arange(1.0, double(n));
        Arr An = 2.0 * xt::cos(xt::linalg::outer(w, cols));
        cached_A = xt::concatenate(xt::xtuple(xt::ones<double>({m, 1}), An), 1);
        cached_n = n;
    }
    const auto& A = cached_A;
    Arr R = xt::linalg::dot(A, r);  // NOQA

    Arr alpha = 0.5 * xt::log(xt::abs(R));

    // find the Hilbert transform
    auto alphatmp = xt::fftw::rfft(alpha);
    // alphatmp(floor(m/2)+1: m) = -alphatmp(floor(m/2)+1: m)
    auto ind = size_t(m) / 2;

    //??? alphatmp[ind:m] = -alphatmp[ind:m];
    xt::view(alphatmp, xt::range(ind, m)) = -xt::view(alphatmp, xt::range(ind, m));
    alphatmp[0] = 0.;
    alphatmp[ind] = 0.;

    // multiply by i*k
    const std::complex<double> j_{0, 1};
    // auto k = xt::fftw::rfftscale<double>(sin.shape()[0], dx);
    // xt::xarray<std::complex<double>> temp= xt::eval(i * alphatmp);
    auto phi = xt::fftw::irfft(xt::xarray<std::complex<double>>(xt::eval(j_ * alphatmp)));

    // now retrieve the original sampling
    // index = find(np.reminder([0:m-1], mult_factor) == 0)
    // auto index = xt::arange(0, m, mult_factor);
    auto alpha1 = xt::view(alpha, xt::range(0, m, mult_factor));
    auto phi1 = xt::view(phi, xt::range(0, m, mult_factor));

    // compute the impulse response (inverse Fourier transform)
    auto h_tmp = xt::exp(alpha1 + j_ * phi1);
    Arr h = xt::fftw::irfft(xt::xarray<std::complex<double>>(xt::eval(h_tmp)));
    return h;
}

/**
 * The `inverse_spectral_fact` function takes an impulse response `h` as input
 * and computes the auto-correlation coefficients `r` that correspond to the
 * given impulse response. It returns the auto-correlation coefficients `r` as
 * an `Arr` object.
 *
 * @param[in] h The parameter `h` is the impulse response, which is a one-dimensional array or
 * vector representing the response of a system to an impulse input.
 *
 * @return The function `inverse_spectral_fact` returns the auto-correlation coefficients `r` as an
 * `Arr` object.
 */
auto inverse_spectral_fact(const Arr& h) -> Arr {
    auto n = h.shape()[0];
    // Use FFT-based autocorrelation (O(n log n) instead of O(n^2))
    // Zero-pad to 2n to avoid circular convolution artifacts
    Arr padded = xt::zeros<double>({2 * n});
    xt::view(padded, xt::range(0, n)) = h;
    auto H = xt::fftw::rfft(padded);
    Arr R = xt::eval(xt::abs(H) * xt::abs(H));            // power spectrum |H|^2
    auto autocorr = xt::fftw::irfft(
        xt::xarray<std::complex<double>>(xt::eval(xt::cast<std::complex<double>>(R))));
    return xt::eval(xt::view(autocorr, xt::range(0, n)) * (2.0 * n));
}
