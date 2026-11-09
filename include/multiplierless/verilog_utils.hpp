#pragma once

/** @file verilog_utils.hpp
 *  @brief Fixed-point and Verilog post-processing helpers.
 */

#include <string>
#include <vector>

/// @brief Coefficients aligned to one shared fixed-point scale.
struct FixedCoefficients {
    std::vector<std::string> digits;  ///< Equal-length integer digit strings
    int max_power = -1;               ///< Exponent of the most significant digit
    int scale_power = 0;              ///< Common fractional width (value * 2^scale_power)
};

/// @brief Right-align fractional CSD digits so every tap shares one scale.
///
/// Each returned string represents `round(value * 2**scale_power)`, where
/// `scale_power` is the common fractional width.
/// @param[in] csd_list CSD strings (with or without a decimal point).
/// @return The aligned digit strings and their common `max_power`.
/// @throws std::invalid_argument If the list is empty or a CSD string is empty.
auto csd_coeffs_to_fixed(const std::vector<std::string>& csd_list) -> FixedCoefficients;

/// @brief Repair csdigit's direct-form module so it compiles.
///
/// csdigit emits a trailing comma after the last port and redeclares each
/// output as `wire signed ... hN = ...`, which conflicts with the port
/// declaration.  This removes the trailing comma and turns those output
/// redeclarations into continuous assignments.
/// @param[in] verilog Module text produced by csd::generate_csd_multipliers().
/// @param[in] port_names Output port names (`h0`, `h1`, ...).
/// @return The repaired module text.
auto fix_direct_verilog(std::string verilog, const std::vector<std::string>& port_names)
    -> std::string;
