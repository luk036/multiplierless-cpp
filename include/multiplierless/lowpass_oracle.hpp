// -*- coding: utf-8 -*-
#pragma once

// #include <limits>
#include <valarray>
#include <xtensor/xarray.hpp>

// static const auto M_PI = std::acos(-1);

// Modified from CVX code by Almir Mutapcic in 2006.
// Adapted in 2010 for impulse response peak-minimization by convex iteration by
// Christine Law.
//
// "FIR Filter Design via Spectral Factorization and Convex Optimization"
// by S.-P. Wu, S. Boyd, and L. Vandenberghe
//
// Designs an FIR lowpass filter using spectral factorization method with
// constraint on maximum passband ripple and stopband attenuation:
//
//   minimize   max |H(w)|                      for w in stopband
//       s.t.   1/delta <= |H(w)| <= delta      for w in passband
//
// We change variables via spectral factorization method and get:
//
//   minimize   max R(w)                          for w in stopband
//       s.t.   (1/delta)**2 <= R(w) <= delta**2  for w in passband
//              R(w) >= 0                         for all w
//
// where R(w) is squared magnitude frequency response
// (and Fourier transform of autocorrelation coefficients r).
// Variables are coeffients r and G = hh' where h is impulse response.
// delta is allowed passband ripple.
// This is a convex problem (can be formulated as an SDP after sampling).

// rand('twister',sum(100*clock))
// randn('state',sum(100*clock))

// *********************************************************************
// filter specs (for a low-pass filter)
// *********************************************************************
// number of FIR coefficients (including zeroth)
struct filter_design_construct {
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

    int N;
    Arr Ap;
    Arr As;
    Arr Anr;
    double Lpsq;
    double Upsq;
    double Spsq;

    filter_design_construct(int argN = 32);
};

// from itertools import chain

/*!
 * @brief Oracle for FIR lowpass filter design.
 *
 *    This example is taken from Almir Mutapcic in 2006:
 *
 *        min   \gamma
 *        s.t.  L^2(\omega) \le R(\omega) \le U^2(\omega), \forall \omega \in
 * [0, \pi] R(\omega) > 0, \forall \omega \in [0, \pi]
 */
class LowpassOracle {
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Vec = std::valarray<double>;
    using ParallelCut = std::pair<Arr, Vec>;

  private:
    size_t _i_Anr{};
    size_t _i_As{};
    size_t _i_Ap{};
    // mutable unsigned int _count{};

    filter_design_construct _Fdc;
    // const Arr& _Ap;
    // const Arr& _As;
    // const Arr& _Anr;
    // double _Lpsq;
    // double _Upsq;

  public:
    bool retry{false}; // ???
    bool more_alt{true};

    /*!
     * @brief Construct a new lowpass oracle object
     *
     * @param[in] Ap
     * @param[in] As
     * @param[in] Anr
     * @param[in] Lpsq
     * @param[in] Upsq
     */
    LowpassOracle(filter_design_construct &&Fdc) : _Fdc{std::move(Fdc)} {}

    /*!
     * @brief
     *
     * @param[in] x
     * @param[in] Spsq
     * @return auto
     */
    auto assess_optim(const Arr &x, double &Spsq)
        -> std::tuple<ParallelCut, bool>;

    /*!
     * @brief
     *
     * @param[in] x
     * @param[in] Spsq
     * @return auto
     */
    auto operator()(const Arr &x, double &Spsq)
        -> std::tuple<ParallelCut, bool> {
        return this->assess_optim(x, Spsq);
    }
};
