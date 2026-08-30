#pragma once

/** @file spectral_fact.hpp
 *  @brief Spectral factorization strategies for minimum-phase FIR filter design.
 */

#include <ellalgo/arr.hpp>

/// @brief Spectral factorization method selector (Strategy pattern).
enum class spectral_method { fft, root };

/// @brief FFT-based spectral factorization (default method).
auto spectral_fact_fft(const Arr& r) -> Arr;

/// @brief Root-finding based spectral factorization (via ginger Aberth solver).
/// @param[in] r         Autocorrelation sequence (top half).
/// @param[in] tolerance Convergence tolerance for the root solver.
/// @return Minimum-phase impulse response coefficients h.
auto spectral_fact_root(const Arr& r, double tolerance) -> Arr;

/// @brief Inverse spectral factorization: autocorrelation of an impulse response.
/// @param[in] h Impulse response coefficients.
/// @return Autocorrelation sequence r (same length as h).
auto inverse_spectral_fact(const Arr& h) -> Arr;

/// @brief Spectral factorization with explicit method selection.
/// @param[in] r         Autocorrelation sequence (top half).
/// @param[in] method    Factorization method (fft or root).
/// @param[in] tolerance Root-solver tolerance (ignored for the fft method).
/// @return Minimum-phase impulse response coefficients h.
inline auto spectral_fact(const Arr& r, spectral_method method, double tolerance) -> Arr {
    switch (method) {
        case spectral_method::root:
            return spectral_fact_root(r, tolerance);
        case spectral_method::fft:
        default:
            return spectral_fact_fft(r);
    }
}

/// @brief Convenience wrapper: default spectral factorization (fft method).
/// @param[in] r Autocorrelation sequence (top half).
/// @return Minimum-phase impulse response coefficients h.
inline auto spectral_fact(const Arr& r) -> Arr { return spectral_fact_fft(r); }
