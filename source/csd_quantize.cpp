/// @file csd_quantize.cpp
/// @brief CSD quantization for multiplierless FIR filter design.
///
/// Implements the greedy double→double CSD quantization declared in
/// <multiplierless/csd_quantize.hpp>. The string-based CSD <-> decimal
/// conversions live in the csd-cpp library (namespace csd).

#include <algorithm>
#include <cmath>
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
