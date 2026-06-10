/// @file csd.cpp
/// @brief CSD wrapper — delegates to csd-cpp library, adds fast double-returning variant.
///
/// Canonical Signed Digit functions for multiplierless FIR filter design.
/// The core CSD ↔ decimal conversions live in the csd-cpp library
/// (namespace csd). This file provides:
///   - Global-namespace wrappers matching the old API used by tests and internals
///   - csd_quantize() — a direct double→double CSD quantization (no string allocs)

#include <cmath>
#include <csd/csd.hpp>
#include <string>

auto csd_quantize(double num, unsigned int nnz) -> double {
    if (num == 0.0) {
        return 0.0;
    }
    auto result = 0.0;
    auto bit_val
        = std::ldexp(1.0, static_cast<int>(std::ceil(std::log2(std::fabs(num) * 1.5))) - 1);
    while (nnz > 0 && std::fabs(num) > 1e-100) {
        if (std::fabs(1.5 * num) > bit_val) {
            auto sgn = (num > 0) ? 1.0 : -1.0;
            result += sgn * bit_val;
            num -= sgn * bit_val;
            --nnz;
        }
        bit_val *= 0.5;
    }
    return result;
}

// =========================================================================
//  Global-namespace wrappers — backward-compatible API for tests / internals
// =========================================================================

auto to_csd(double num, int places) -> std::string {
    if (num == 0.0) {
        return "0";
    }
    return csd::to_csd(num, places);
}

auto to_decimal(const std::string& csd_str) -> double { return csd::to_decimal(csd_str.c_str()); }

auto to_csdnnz(double num, unsigned int nnz) -> std::string { return csd::to_csdnnz(num, nnz); }
