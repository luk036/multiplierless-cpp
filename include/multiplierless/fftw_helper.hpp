#pragma once

/// @file fftw_helper.hpp
/// @brief FFTW wrappers and complex array ops for spectral factorization.
///
/// Provides FFT/IFFT wrappers around the FFTW3 library plus complex array
/// arithmetic operators. These are kept separate from the core Arr type
/// (which now lives in ellalgo-cpp/include/ellalgo/arr.hpp).

#include <fftw3.h>

#include <complex>
#include <ellalgo/arr.hpp>
#include <vector>

/** @brief Complex array type used for FFT operations. */
using CArr = std::vector<std::complex<double>>;

/// @brief Cast a real Arr to a complex CArr.
/// @param a Real-valued input array.
/// @return Complex array with zero imaginary part.
inline CArr cast_to_complex(const Arr& a) {
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = std::complex<double>(a(i), 0.0);
    return out;
}

/// @brief Real-input FFT (r2c 1D).
/// Computes the one-dimensional FFT of a real-valued array using FFTW's r2c transform.
/// @param in Real-valued input array.
/// @return Complex frequency-domain representation (n/2 + 1 elements).
inline CArr rfft(const Arr& in) {
    size_t n = in.size();
    size_t n_out = n / 2 + 1;
    CArr out(n_out);
    auto* plan = fftw_plan_dft_r2c_1d(static_cast<int>(n), const_cast<double*>(in.data()),
                                      reinterpret_cast<fftw_complex*>(out.data()), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return out;
}

/// @brief Inverse real FFT (c2r 1D).
/// Computes the inverse FFT of a complex array back to a real-valued result,
/// with a 1/N scaling applied to match the MATLAB/numpy convention.
/// @param in Complex frequency-domain input (n/2 + 1 elements from rfft).
/// @param expected_size Expected output size (original N before rfft).
/// @return Real-valued time-domain array of length expected_size.
inline Arr irfft(const CArr& in, size_t expected_size) {
    Arr out(expected_size);
    auto* plan = fftw_plan_dft_c2r_1d(
        static_cast<int>(expected_size),
        reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(in.data())), out.data(),
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (size_t i = 0; i < expected_size; ++i) out(i) /= static_cast<double>(expected_size);
    return out;
}

/// @brief Full complex-to-complex FFT (like numpy.fft.fft).
/// No scaling is applied to the output.
/// @param in Complex-valued input array.
/// @return Complex frequency-domain representation (same size as input).
inline CArr fft(const CArr& in) {
    size_t n = in.size();
    CArr out(n);
    auto* plan = fftw_plan_dft_1d(
        static_cast<int>(n),
        reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(in.data())),
        reinterpret_cast<fftw_complex*>(out.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return out;
}

/// @brief Full complex-to-real IFFT (like numpy.fft.ifft, taking real part).
/// The result is scaled by 1/N to match the numpy convention.
/// @param in Complex-valued frequency-domain input.
/// @return Real-valued time-domain array (same size as input).
inline Arr ifft(const CArr& in) {
    size_t n = in.size();
    CArr tmp(n);
    auto* plan = fftw_plan_dft_1d(
        static_cast<int>(n),
        reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(in.data())),
        reinterpret_cast<fftw_complex*>(tmp.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    Arr out(n);
    for (size_t i = 0; i < n; ++i) out(i) = tmp[i].real() / static_cast<double>(n);
    return out;
}

/// @brief Extract the real part of a complex array.
/// @param a Complex-valued input array.
/// @return Real-valued array containing the real components.
inline Arr real(const CArr& a) {
    Arr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out(i) = a[i].real();
    return out;
}

/// @defgroup complex_ops Complex array arithmetic
/// @brief Element-wise arithmetic operators for CArr (complex vector).
/// @{
/// @brief Scalar multiplication of complex array.
inline CArr operator*(const std::complex<double>& s, const CArr& a) {
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = s * a[i];
    return out;
}

/// @brief Element-wise complex + complex addition.
inline CArr operator+(const CArr& a, const CArr& b) {
    assert(a.size() == b.size());
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = a[i] + b[i];
    return out;
}

/// @brief Element-wise real + complex addition.
inline CArr operator+(const Arr& a, const CArr& b) {
    assert(a.size() == b.size());
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = std::complex<double>(a(i), 0.0) + b[i];
    return out;
}

/// @brief Element-wise complex exponential.
inline CArr exp(const CArr& a) {
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = std::exp(a[i]);
    return out;
}

/// @brief Element-wise complex magnitude.
inline Arr abs(const CArr& a) {
    Arr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out(i) = std::abs(a[i]);
    return out;
}
/// @}
