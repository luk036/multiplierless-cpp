// -*- coding: utf-8 -*-
#pragma once

// Disable svector on macOS to avoid Clang template ambiguity issues
// where long and unsigned long are both 64-bit
#ifdef __APPLE__
#define XTENSOR_DISABLE_SVECTOR 1
#endif

// #include <limits>
#include <valarray>
#include <xtensor/xarray.hpp>

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
/// Filter design construct containing all parameters for lowpass filter design.
/// This struct holds the filter order, passband/stopband specifications, and
/// squared constraints used in the spectral factorization method for FIR filter design.
/// @param argN The filter order (number of FIR coefficients including zeroth)
struct filter_design_construct {
    using Arr = xt::xarray<double>;

    int N;                      ///< Filter order (number of FIR coefficients including zeroth)
    Arr Ap;                     ///< Passband constraint matrix
    Arr As;                     ///< Stopband constraint matrix
    Arr Anr;                    ///< Non-redundant constraint matrix
    double Lpsq;                ///< Lower bound squared for passband (1/delta^2)
    double Upsq;                ///< Upper bound squared for passband (delta^2)
    double Spsq;                ///< Stopband attenuation squared

    /// Construct a filter_design_construct with default parameters for the given filter order.
    /// @param argN The filter order (N+1 coefficients will be generated)
    explicit filter_design_construct(int argN = 32);
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
    using Arr = xt::xarray<double>;
    using Vec = std::valarray<double>;
    using ParallelCut = std::pair<Arr, Vec>;

    size_t _i_Anr{};
    size_t _i_As{};
    size_t _i_Ap{};

    filter_design_construct _Fdc;

  public:
    // bool retry{false};  // ???
    // bool more_alt{true};

    /*!
     * @brief Construct a new lowpass oracle object
     *
     * @param[in] Fdc A filter_design_construct object containing all necessary parameters
     *                for the lowpass filter design (Ap, As, Anr, Lpsq, Upsq, etc.)
     */
    explicit LowpassOracle(filter_design_construct &&Fdc) : _Fdc{std::move(Fdc)} {}

    /*!
     * @brief Assess the optimization problem for the given filter coefficients.
     *
     * This function evaluates the constraints and computes the cutting plane
     * for the FIR lowpass filter design optimization problem. It checks the
     * non-negative real constraint, passband constraints, stopband constraints,
     * and returns the appropriate gradient and objective function values.
     *
     * @param[in] x The filter coefficients (autocorrelation coefficients r).
     * @param[in,out] Spsq On input, the initial stopband attenuation squared target.
     *                      On output, the achieved maximum stopband value.
     *
     * @return A tuple containing:
     *         - ParallelCut: Pair of gradient and objective function values
     *         - bool: True if optimal solution found, false if more iterations needed
     */
    auto assess_optim(const Arr &x, double &Spsq) -> std::tuple<ParallelCut, bool>;

    /*!
     * @brief Operator function for optimization assessment.
     *
     * This is a convenience function that forwards to assess_optim.
     *
     * @param[in] x The filter coefficients (autocorrelation coefficients r).
     * @param[in,out] Spsq On input, the initial stopband attenuation squared target.
     *                      On output, the achieved maximum stopband value.
     *
     * @return A tuple containing:
     *         - ParallelCut: Pair of gradient and objective function values
     *         - bool: True if optimal solution found, false if more iterations needed
     */
    auto operator()(const Arr &x, double &Spsq) -> std::tuple<ParallelCut, bool> {
        return this->assess_optim(x, Spsq);
    }
};
