#pragma once

/// FFTW wrappers and complex array ops, kept separate from the core Arr
/// (which now lives in ellalgo-cpp/include/ellalgo/arr.hpp).

#include <complex>
#include <fftw3.h>
#include <vector>

#include <ellalgo/arr.hpp>

// Complex array type
using CArr = std::vector<std::complex<double>>;

// Real → complex cast
inline CArr cast_to_complex(const Arr& a) {
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = std::complex<double>(a(i), 0.0);
    return out;
}

// Real → complex FFT
inline CArr rfft(const Arr& in) {
    size_t n = in.size();
    size_t n_out = n / 2 + 1;
    CArr out(n_out);
    auto* plan = fftw_plan_dft_r2c_1d(
        static_cast<int>(n),
        const_cast<double*>(in.data()),
        reinterpret_cast<fftw_complex*>(out.data()),
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return out;
}

// Complex → real FFT
inline Arr irfft(const CArr& in, size_t expected_size) {
    Arr out(expected_size);
    auto* plan = fftw_plan_dft_c2r_1d(
        static_cast<int>(expected_size),
        reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(in.data())),
        out.data(),
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for (size_t i = 0; i < expected_size; ++i) out(i) /= static_cast<double>(expected_size);
    return out;
}

// Full complex → complex FFT (like np.fft.fft) — no scaling
inline CArr fft(const CArr& in) {
    size_t n = in.size();
    CArr out(n);
    auto* plan = fftw_plan_dft_1d(
        static_cast<int>(n),
        reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(in.data())),
        reinterpret_cast<fftw_complex*>(out.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return out;
}

// Full complex → real IFFT (like np.fft.ifft, taking real part)
inline Arr ifft(const CArr& in) {
    size_t n = in.size();
    CArr tmp(n);
    auto* plan = fftw_plan_dft_1d(
        static_cast<int>(n),
        reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(in.data())),
        reinterpret_cast<fftw_complex*>(tmp.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    Arr out(n);
    for (size_t i = 0; i < n; ++i) out(i) = tmp[i].real() / static_cast<double>(n);
    return out;
}

// Real part of complex array
inline Arr real(const CArr& a) {
    Arr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out(i) = a[i].real();
    return out;
}

// Complex element-wise operations
inline CArr operator*(const std::complex<double>& s, const CArr& a) {
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = s * a[i];
    return out;
}

inline CArr operator+(const CArr& a, const CArr& b) {
    assert(a.size() == b.size());
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = a[i] + b[i];
    return out;
}

inline CArr operator+(const Arr& a, const CArr& b) {
    assert(a.size() == b.size());
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = std::complex<double>(a(i), 0.0) + b[i];
    return out;
}

inline CArr exp(const CArr& a) {
    CArr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out[i] = std::exp(a[i]);
    return out;
}

inline Arr abs(const CArr& a) {
    Arr out(a.size());
    for (size_t i = 0; i < a.size(); ++i) out(i) = std::abs(a[i]);
    return out;
}
