#pragma once

/** @file csd_quantize.hpp
 *  @brief CSD quantization helpers shared by the optimizer and the CLI.
 */

#include <string>
#include <utility>
#include <vector>

/// @brief A quantized value together with the CSD digits that produced it.
///
/// Each digit is a pair `(exponent, sign)` where `sign` is `+1` or `-1` and
/// the digit contributes `sign * 2^exponent`.  The digits sum exactly to
/// `value`, so callers can render a CSD string without re-quantizing.
struct CsdParts {
    double value = 0.0;
    std::vector<std::pair<int, int>> digits;
};

/// @brief Greedy CSD quantization returning both the value and its digits.
/// @param[in] num Input value.
/// @param[in] nnz Maximum number of non-zero digits.
/// @return The quantized value and the chosen `(exponent, sign)` digits.
auto csd_quantize_parts(double num, unsigned int nnz) -> CsdParts;

/// @brief Greedy double-to-double CSD quantization.
/// @param[in] num Input value.
/// @param[in] nnz Maximum number of non-zero digits.
/// @return The quantized value.
auto csd_quantize(double num, unsigned int nnz) -> double;

/// @brief Render CSD digits as an exact `integer.fraction` string.
///
/// A decimal point and at least one fractional digit are always emitted, so
/// the output schema holds even for values needing no fraction (0 or ±1).
/// @param[in] digits Digits returned by csd_quantize_parts().
/// @return The CSD string.
auto format_csd(const std::vector<std::pair<int, int>>& digits) -> std::string;

/// @brief Exact numeric value of a CSD string (decimal point optional).
/// @param[in] csd CSD string using '+', '-', '0', and optionally '.'.
/// @return The represented value.
auto csd_to_float(const std::string& csd) -> double;
