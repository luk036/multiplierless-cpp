#include <cmath>   // for cos, M_PI
#include <complex> // for complex, operator*, operator-
#include <tuple>   // for tuple

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#include <xtensor-blas/xlinalg.hpp> // for dot
#include <xtensor-fftw/basic.hpp>   // for irfft, rfft
#include <xtensor/xaccessible.hpp>  // for xaccessible
#include <xtensor/xarray.hpp>       // for xarray_container
#include <xtensor/xbroadcast.hpp>   // for xbroadcast
#include <xtensor/xbuilder.hpp>     // for zeros, concatenate, linspace
#include <xtensor/xcontainer.hpp>   // for xcontainer, xcontainer<>::...
#include <xtensor/xeval.hpp>        // for eval
#include <xtensor/xexception.hpp>   // for throw_concatenate_error
#include <xtensor/xfunction.hpp>    // for xfunction
#include <xtensor/xgenerator.hpp>   // for xgenerator
#include <xtensor/xiterator.hpp>    // for linear_begin
#include <xtensor/xlayout.hpp>      // for layout_type, layout_type::...
#include <xtensor/xmath.hpp>        // for abs, exp, log, sum, abs_fun
#include <xtensor/xoperation.hpp>   // for xfunction_type_t, operator*
#include <xtensor/xreducer.hpp>     // for xreducer
#include <xtensor/xslice.hpp>       // for range, xtuph, _
#include <xtensor/xutils.hpp>       // for accumulate
#include <xtensor/xview.hpp>        // for xview, view

using Arr = xt::xarray<double, xt::layout_type::row_major>;

/* The `spectral_fact` function is performing spectral factorization using the
Kolmogorov 1939 approach. It takes an input vector `r` which represents the
top-half of the auto-correlation coefficients. It computes the minimum-phase
impulse response `h` that satisfies the given auto-correlation. The function
returns the impulse response `h` as an `Arr` object. */

/**
 * @brief Spectral factorization
 *
 *    Spectral factorization using Kolmogorov 1939 approach.
 *      (code follows pp. 232-233, Signal Analysis, by A. Papoulis)
 *
 *      Computes the minimum-phase impulse response which satisfies
 *      given auto-correlation.
 *
 *      Input:
 *        r: top-half of the auto-correlation coefficients
 *           starts from 0th element to end of the auto-corelation
 *           should be passed in as a column vector
 *      Output
 *        h: impulse response that gives the desired auto-correlation
 *
 * @return auto
 */
auto spectral_fact(const Arr &r) -> Arr {
  // length of the impulse response sequence
  const auto n = int(r.shape()[0]);

  // over-sampling factor
  const auto mult_factor = 100; // should have mult_factor*(n) >> n
  const auto m = mult_factor * n;
  // const auto PI = std::acos(-1);

  // computation method:
  // H(exp(jTw)) = alpha(w) + j*phi(w)
  // where alpha(w) = 1/2*ln(R(w)) and phi(w) = Hilbert_trans(alpha(w))

  // compute 1/2*ln(R(w))
  // w = 2*pi*[0:m-1]/m
  Arr w = xt::linspace<double>(0, 2 * M_PI, size_t(m));

  // R = [ones(m, 1) 2*cos(kron(w', [1:n-1]))]*r
  // Arr Bn = xt::linalg::outer(w, xt::arange(1, n));
  // Arr An = 2 * xt::cos(Bn);
  // // Arr R = np.hstack((np.ones((m, 1)), An)).dot(r)
  // Arr A = xt::concatenate(xt::xtuple(xt::ones<double>({m, 1}), An), 1);
  auto An = Arr(xt::zeros<double>({m, n - 1}));
  for (auto i = 0; i != m; ++i) {
    for (auto j = 0; j != n - 1; ++j) {
      An(i, j) = 2.0 * std::cos(w(i) * (j + 1));
    }
  }
  Arr A = xt::concatenate(xt::xtuple(xt::ones<double>({m, 1}), An), 1);
  Arr R = xt::linalg::dot(A, r); // NOQA

  Arr alpha = 0.5 * xt::log(xt::abs(R));

  // find the Hilbert transform
  auto alphatmp = xt::fftw::rfft(alpha);
  // alphatmp(floor(m/2)+1: m) = -alphatmp(floor(m/2)+1: m)
  auto ind = size_t(m) / 2;

  //??? alphatmp[ind:m] = -alphatmp[ind:m];
  xt::view(alphatmp, xt::range(ind, m)) =
      -xt::view(alphatmp, xt::range(ind, m));
  alphatmp[0] = 0.;
  alphatmp[ind] = 0.;

  // multiply by i*k
  const std::complex<double> j_{0, 1};
  // auto k = xt::fftw::rfftscale<double>(sin.shape()[0], dx);
  // xt::xarray<std::complex<double>> temp= xt::eval(i * alphatmp);
  auto phi = xt::fftw::irfft(xt::eval(j_ * alphatmp));

  // now retrieve the original sampling
  // index = find(np.reminder([0:m-1], mult_factor) == 0)
  // auto index = xt::arange(0, m, mult_factor);
  auto alpha1 = xt::view(alpha, xt::range(0, m, mult_factor));
  auto phi1 = xt::view(phi, xt::range(0, m, mult_factor));

  // compute the impulse response (inverse Fourier transform)
  Arr h = xt::fftw::irfft(xt::eval(xt::exp(alpha1 + j_ * phi1)));
  return h;
}

/**
 * @brief
 *
 * The `inverse_spectral_fact` function takes an impulse response `h` as input and
 * computes the auto-correlation coefficients `r` that correspond to the given
 * impulse response. It returns the auto-correlation coefficients `r` as an `Arr`
 * object.
 *
 * @param h
 * @return Arr
 */
auto inverse_spectral_fact(const Arr &h) -> Arr {
  auto n = h.shape()[0];
  auto r = Arr{xt::zeros<double>({n})};
  using xt::placeholders::_;
  for (auto t = 0U; t != n; ++t) {
    r(t) = xt::sum(xt::view(h, xt::range(t, _)) *
                   xt::view(h, xt::range(_, n - t)))();
  }
  return r;
}

// if __name__ == "__main__":
//     r = np.random.rand(20)
//     h = spectral_fact(r)
//     print(h)
