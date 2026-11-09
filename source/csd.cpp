/// @file csd.cpp
/// @brief CSD wrapper — delegates to csd-cpp library, adds fast double-returning variant.
///
/// Canonical Signed Digit functions for multiplierless FIR filter design.
/// The core CSD ↔ decimal conversions live in the csd-cpp library
/// (namespace csd). This file provides:
///   - Global-namespace wrappers matching the old API used by tests and internals
///   - csd_quantize() — a direct double→double CSD quantization (no string allocs)

#include <algorithm>
#include <cmath>
#include <csd/csd.hpp>
#include <map>
#include <multiplierless/csd_quantize.hpp>
#include <string>

auto csd_quantize_parts(double num, unsigned int nnz) -> CsdParts {
    if (num == 0.0 || nnz == 0) {
        return {};
    }
    auto result = 0.0;
    std::vector<std::pair<int, int>> digits;
    auto exp = static_cast<int>(std::ceil(std::log2(std::fabs(num) * 1.5))) - 1;
    while (nnz > 0 && std::fabs(num) > 1e-100) {
        auto const bit_val = std::ldexp(1.0, exp);
        if (std::fabs(1.5 * num) > bit_val) {
            auto const sgn = (num > 0) ? 1 : -1;
            result += sgn * bit_val;
            digits.emplace_back(exp, sgn);
            num -= sgn * bit_val;
            --nnz;
        }
        --exp;
    }
    return {result, std::move(digits)};
}

auto csd_quantize(double num, unsigned int nnz) -> double {
    return csd_quantize_parts(num, nnz).value;
}

auto format_csd(const std::vector<std::pair<int, int>>& digits) -> std::string {
    if (digits.empty()) {
        return "0.";
    }
    std::map<int, char> by_exp;
    auto max_exp = 0;
    auto min_exp = -1;
    for (auto const& [exp, sign] : digits) {
        by_exp[exp] = (sign > 0) ? '+' : '-';
        max_exp = std::max(max_exp, exp);
        min_exp = std::min(min_exp, exp);
    }
    std::string out;
    for (auto exp = max_exp; exp >= 0; --exp) {
        auto const it = by_exp.find(exp);
        out += (it == by_exp.end()) ? '0' : it->second;
    }
    out += '.';
    for (auto exp = -1; exp >= min_exp; --exp) {
        auto const it = by_exp.find(exp);
        out += (it == by_exp.end()) ? '0' : it->second;
    }
    return out;
}

auto csd_to_float(const std::string& csd) -> double {
    auto const dot = csd.find('.');
    auto const int_part = (dot == std::string::npos) ? csd : csd.substr(0, dot);
    auto const frac_part = (dot == std::string::npos) ? std::string{} : csd.substr(dot + 1);
    auto value = 0.0;
    for (size_t i = 0; i < int_part.size(); ++i) {
        auto const exp = static_cast<int>(int_part.size() - 1 - i);
        if (int_part[i] == '+') value += std::ldexp(1.0, exp);
        if (int_part[i] == '-') value -= std::ldexp(1.0, exp);
    }
    for (size_t i = 0; i < frac_part.size(); ++i) {
        auto const exp = -static_cast<int>(i + 1);
        if (frac_part[i] == '+') value += std::ldexp(1.0, exp);
        if (frac_part[i] == '-') value -= std::ldexp(1.0, exp);
    }
    return value;
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
