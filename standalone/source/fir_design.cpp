#include <algorithm>
#include <cmath>
#include <csd/csd_multiplier.hpp>
#include <ellalgo/arr.hpp>
#include <ellalgo/cutting_plane.hpp>
#include <ellalgo/ell.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <multiplierless/lowpass_oracle.hpp>
#include <multiplierless/lowpass_oracle_q.hpp>
#include <nlohmann/json.hpp>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using json = nlohmann::json;

// ============================================================
//  Transpose-form FIR filter Verilog generator (self-contained)
// ============================================================
namespace {

    // Count non-zero CSD digits ('+' or '-') in a string
    auto count_nnz(const std::string& s) -> int {
        int n = 0;
        for (auto c : s) {
            if (c == '+' || c == '-') ++n;
        }
        return n;
    }

    // Build a flat Verilog expression for a CSD range using x_shift references
    auto build_range_expr(const std::string& csd_str, size_t start, size_t length, int max_power)
        -> std::string {
        std::string expr;
        bool first = true;
        for (size_t i = start; i < start + length && i < csd_str.size(); ++i) {
            auto const power = max_power - static_cast<int>(i);
            switch (csd_str[i]) {
                case '+':
                    if (first) {
                        expr += "x_shift" + std::to_string(power);
                        first = false;
                    } else {
                        expr += " + x_shift" + std::to_string(power);
                    }
                    break;
                case '-':
                    if (first) {
                        expr += "-x_shift" + std::to_string(power);
                        first = false;
                    } else {
                        expr += " - x_shift" + std::to_string(power);
                    }
                    break;
                default:
                    break;
            }
        }
        return expr;
    }

    // Find non-overlapping occurrences of |pattern| in |csd_str|
    auto find_pattern_occurrences(const std::string& csd_str, const std::string& pattern)
        -> std::vector<size_t> {
        std::vector<size_t> positions;
        size_t pos = 0;
        while ((pos = csd_str.find(pattern, pos)) != std::string::npos) {
            positions.push_back(pos);
            pos += pattern.size();
        }
        return positions;
    }

    // Find substrings (NNZ >= 2) shared by >= 2 different CSD strings
    auto find_cross_patterns(const std::vector<std::string>& csd_list)
        -> std::map<std::string, std::vector<std::pair<int, int>>> {
        std::map<std::string, std::vector<std::pair<int, int>>> patterns;
        for (int ci = 0; ci < static_cast<int>(csd_list.size()); ++ci) {
            auto const& csd = csd_list[ci];
            auto const n = static_cast<int>(csd.size());
            for (int i = 0; i < n; ++i) {
                for (int j = i + 2; j <= n; ++j) {
                    auto sub = csd.substr(i, j - i);
                    if (count_nnz(sub) >= 2) {
                        patterns[sub].emplace_back(ci, i);
                    }
                }
            }
        }
        // Keep only patterns crossing >= 2 different CSD strings
        auto it = patterns.begin();
        while (it != patterns.end()) {
            std::set<int> unique_idx;
            for (auto const& occ_pair : it->second) {
                unique_idx.insert(occ_pair.first);
            }
            if (unique_idx.size() < 2) {
                it = patterns.erase(it);
            } else {
                ++it;
            }
        }
        return patterns;
    }

    // Build coefficient expression using CSE wire + flat gap terms
    auto build_coeff_expr(const std::string& csd, int max_power, const std::string& pattern,
                          int base_pos, const std::string& cse_name) -> std::string {
        if (pattern.empty()) {
            return build_range_expr(csd, 0, csd.size(), max_power);
        }
        auto const positions = find_pattern_occurrences(csd, pattern);
        std::vector<std::string> parts;
        size_t cur = 0;
        for (auto const pos : positions) {
            if (pos > cur) {
                auto gap = build_range_expr(csd, cur, pos - cur, max_power);
                if (!gap.empty()) parts.push_back(gap);
            }
            auto const shift = static_cast<int>(pos) - base_pos;
            if (shift == 0) {
                parts.push_back(cse_name);
            } else {
                parts.push_back("(" + cse_name + " >>> " + std::to_string(shift) + ")");
            }
            cur = pos + pattern.size();
        }
        if (cur < csd.size()) {
            auto gap = build_range_expr(csd, cur, csd.size() - cur, max_power);
            if (!gap.empty()) parts.push_back(gap);
        }
        if (parts.empty()) return {};
        std::string result = parts[0];
        for (size_t i = 1; i < parts.size(); ++i) {
            result += " + " + parts[i];
        }
        return result;
    }

    // Generate transpose-form FIR filter Verilog with cross-CSE
    auto generate_transpose_form_verilog(const std::vector<csd::MultiplierSpec>& coeffs,
                                         const std::string& module_name) -> std::string {
        if (coeffs.empty()) {
            throw std::invalid_argument("At least one coefficient is required");
        }
        auto const input_width = coeffs[0].input_width;
        auto const max_power = coeffs[0].max_power;
        auto const N = static_cast<int>(coeffs.size());
        auto const output_width = input_width + max_power;

        // Validation
        for (auto const& spec : coeffs) {
            if (spec.input_width != input_width || spec.max_power != max_power) {
                throw std::invalid_argument(
                    "All coefficients must share input_width and max_power");
            }
            auto const len = static_cast<int>(spec.csd.size());
            if (len != max_power + 1) {
                throw std::invalid_argument("CSD length mismatch for '" + spec.name + "'");
            }
            for (auto c : spec.csd) {
                if (c != '+' && c != '-' && c != '0') {
                    throw std::invalid_argument("CSD string can only contain '+', '-', or '0'");
                }
            }
        }

        // Collect shift powers
        std::set<int, std::greater<int>> all_powers;
        for (auto const& spec : coeffs) {
            for (int i = 0; i < static_cast<int>(spec.csd.size()); ++i) {
                if (spec.csd[i] != '0') {
                    all_powers.insert(max_power - i);
                }
            }
        }

        // Cross-CSE pattern detection
        std::vector<std::string> csd_strings;
        csd_strings.reserve(N);
        for (auto const& spec : coeffs) {
            csd_strings.push_back(spec.csd);
        }
        auto const cross = find_cross_patterns(csd_strings);

        std::string best_pattern;
        std::vector<std::pair<int, int>> best_occurrences;
        int best_score = 0;
        for (auto const& kv : cross) {
            auto const& pat = kv.first;
            auto const& occ = kv.second;
            auto const nnz = count_nnz(pat);
            auto const score_val = (nnz - 1) * (static_cast<int>(occ.size()) - 1);
            if (score_val > best_score) {
                best_score = score_val;
                best_pattern = pat;
                best_occurrences = occ;
            }
        }

        int cse_base_pos = 0;
        if (!best_pattern.empty()) {
            cse_base_pos = best_occurrences[0].second;
            for (auto const& occ_pair : best_occurrences) {
                if (occ_pair.second < cse_base_pos) {
                    cse_base_pos = occ_pair.second;
                }
            }
        }
        std::set<int> cse_coeffs;
        for (auto const& occ_pair : best_occurrences) {
            cse_coeffs.insert(occ_pair.first);
        }

        // Build Verilog module
        std::string verilog;
        verilog += "\nmodule " + module_name + " (";
        verilog += "\n    input clk,";
        verilog += "\n    input rst_n,";
        verilog += "\n    input signed [" + std::to_string(input_width - 1) + ":0] x,";
        verilog += "\n    output signed [" + std::to_string(output_width - 1) + ":0] y";
        verilog += "\n);";

        if (!all_powers.empty()) {
            verilog += "\n\n    // Shifted versions of input";
            for (auto p : all_powers) {
                verilog += "\n    wire signed [" + std::to_string(output_width - 1) + ":0] x_shift"
                           + std::to_string(p) + " = x <<< " + std::to_string(p) + ";";
            }
        }

        if (!best_pattern.empty()) {
            auto const cse_expr
                = build_range_expr(best_pattern, 0, best_pattern.size(), max_power - cse_base_pos);
            verilog += "\n\n    // Cross-CSE: shared pattern \"" + best_pattern + "\"";
            verilog += "\n    wire signed [" + std::to_string(output_width - 1)
                       + ":0] _cse_0 = " + cse_expr + ";";
        }

        verilog += "\n\n    // Transpose-form pipeline registers";
        for (int idx = 0; idx < N; ++idx) {
            verilog += "\n    reg signed [" + std::to_string(output_width - 1) + ":0] sum"
                       + std::to_string(idx) + ";";
        }

        verilog += "\n\n    always @(posedge clk or negedge rst_n) begin";
        verilog += "\n        if (!rst_n) begin";
        for (int idx = 0; idx < N; ++idx) {
            verilog += "\n            sum" + std::to_string(idx) + " <= 0;";
        }
        verilog += "\n        end else begin";

        // Coefficients in REVERSE order (canonical transpose form)
        for (int idx = 0; idx < N; ++idx) {
            auto const coeff_idx = N - 1 - idx;
            auto const& spec = coeffs[coeff_idx];
            bool has_cse = !best_pattern.empty() && (cse_coeffs.count(coeff_idx) != 0U);
            auto const expr = has_cse ? build_coeff_expr(spec.csd, max_power, best_pattern,
                                                         cse_base_pos, "_cse_0")
                                      : build_coeff_expr(spec.csd, max_power, {}, 0, {});

            if (idx == 0) {
                if (expr.empty()) {
                    verilog += "\n            sum0 <= 0;";
                } else {
                    verilog += "\n            sum0 <= " + expr + ";";
                }
            } else {
                if (expr.empty()) {
                    verilog += "\n            sum" + std::to_string(idx) + " <= sum"
                               + std::to_string(idx - 1) + ";";
                } else {
                    verilog += "\n            sum" + std::to_string(idx) + " <= sum"
                               + std::to_string(idx - 1) + " + " + expr + ";";
                }
            }
        }

        verilog += "\n        end";
        verilog += "\n    end";

        verilog += "\n\n    assign y = sum" + std::to_string(N - 1) + ";";
        verilog += "\nendmodule\n";
        return verilog;
    }

    // Fix missing commas between port declarations (known csd-cpp issue)
    auto fix_verilog_ports(const std::string& verilog) -> std::string {
        // Collect all lines
        std::vector<std::string> lines;
        size_t pos = 0;
        while (pos < verilog.size()) {
            auto eol = verilog.find('\n', pos);
            lines.push_back(verilog.substr(pos, eol - pos));
            pos = (eol != std::string::npos) ? eol + 1 : verilog.size();
        }

        // Identify port lines (between module header and ");")
        int module_start = -1;
        int paren_end = -1;
        std::vector<int> port_lines;
        for (int i = 0; i < static_cast<int>(lines.size()); ++i) {
            auto& line = lines[i];
            auto trimmed = line;
            trimmed.erase(0, trimmed.find_first_not_of(" \t\r"));
            if (module_start < 0 && trimmed.find("module ") == 0
                && line.find('(') != std::string::npos) {
                module_start = i;
            }
            if (module_start >= 0 && paren_end < 0) {
                if (trimmed == ");") {
                    paren_end = i;
                } else if (i != module_start) {
                    port_lines.push_back(i);
                }
            }
        }

        // Add commas to all port lines except the last
        for (size_t idx = 0; idx < port_lines.size(); ++idx) {
            if (idx == port_lines.size() - 1) continue;  // skip last
            auto& line = lines[port_lines[idx]];
            auto comment_pos = line.find("//");
            auto code_end = (comment_pos != std::string::npos) ? comment_pos : line.size();
            auto code = line.substr(0, code_end);
            while (!code.empty()
                   && (code.back() == ' ' || code.back() == '\t' || code.back() == '\r')) {
                code.pop_back();
            }
            if (!code.empty() && code.back() != ',') {
                line.insert(code_end, ",");
            }
        }

        // Rejoin
        std::string result;
        for (size_t i = 0; i < lines.size(); ++i) {
            result += lines[i];
            if (i + 1 < lines.size()) result += '\n';
        }
        return result;
    }
}  // anonymous namespace

extern auto csd_quantize(double num, unsigned int nnz) -> double;
extern auto spectral_fact_fft(const Arr& r) -> Arr;
extern auto spectral_fact_root(const Arr& r, double tolerance) -> Arr;
extern auto to_csdnnz(double num, unsigned int nnz) -> std::string;

/**
 * @brief FIR filter design tool — main entry point.
 *
 * Reads a JSON filter specification, runs ellipsoid-method optimisation with
 * CSD quantisation, performs spectral factorisation, and outputs the
 * multiplierless FIR coefficients in JSON format (with optional Verilog
 * generation for the CSD multiplier hardware).
 *
 * @param[in] argc Argument count.
 * @param[in] argv Argument vector. Usage: fir_design <filter_spec.json>
 * @return 0 on success, 1 on error.
 */
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filter_spec.json>\n";
        std::cerr << "  Reads filter specifications from a JSON file and\n";
        std::cerr << "  outputs multiplierless FIR coefficients in CSD format.\n";
        return 1;
    }

    std::ifstream ifs(argv[1]);
    if (!ifs) {
        std::cerr << "Error: cannot open " << argv[1] << '\n';
        return 1;
    }

    json spec;
    try {
        ifs >> spec;
    } catch (const json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << '\n';
        return 1;
    }

    auto filter_order = spec.value("filter_order", 32);
    auto wpass_norm = spec.value("passband_edge", 0.12);
    auto wstop_norm = spec.value("stopband_edge", 0.20);
    auto passband_ripple = spec.value("passband_ripple", 0.125);
    auto stopband_attn = spec.value("stopband_attenuation", 0.125);
    auto csd_nnz = spec.value("csd_nnz", 7);
    auto disc_factor = spec.value("discretization_factor", 15);
    auto max_iters = spec.value("max_iters", static_cast<size_t>(50000));
    auto tolerance = spec.value("tolerance", 1e-14);
    auto ell_radius = spec.value("ellipsoid_radius", 40.0);
    auto parallel_cut = spec.value("parallel_cut", true);
    auto gen_verilog = spec.contains("verilog");

    auto Fdc = filter_design_construct(filter_order, wpass_norm, wstop_norm, passband_ripple,
                                       stopband_attn, disc_factor);
    auto Spsq = Fdc.Spsq;

    auto omega = LowpassOracleQ(csd_nnz, LowpassOracle(std::move(Fdc)));

    auto r0 = zeros(filter_order);
    auto ellip = Ell<Arr>(ell_radius, r0);
    ellip.set_use_parallel_cut(parallel_cut);

    Options options(max_iters, tolerance);

    auto [r, num_iters] = cutting_plane_optim_q(omega, ellip, Spsq, options);

    if (r.size() == 0) {
        std::cerr << "Optimization failed — no feasible solution found after " << num_iters
                  << " iterations.\n";
        return 1;
    }

    auto spectral_method = spec.value("spectral_method", std::string("fft"));
    auto root_tol = spec.value("root_tolerance", 1e-8);
    Arr h;
    if (spectral_method == "fft") {
        h = spectral_fact_fft(r);
    } else {
        h = spectral_fact_root(r, root_tol);
    }

    json output;
    output["filter_order"] = filter_order;
    output["csd_nnz"] = csd_nnz;
    output["iterations"] = num_iters;
    output["spectral_method"] = spectral_method;
    output["coefficients"] = json::array();

    for (size_t i = 0; i < h.size(); ++i) {
        auto coeff = csd_quantize(h(i), csd_nnz);
        auto csd_str = to_csdnnz(h(i), csd_nnz);
        output["coefficients"].push_back({{"index", i}, {"value", coeff}, {"csd", csd_str}});
    }

    if (gen_verilog) {
        auto& vl = spec["verilog"];
        auto input_width = vl.value("input_width", 16);
        auto module_name = vl.value("module_name", "fir_filter");
        auto verilog_form = vl.value("form", std::string("transpose"));

        std::vector<std::string> csd_strings;
        for (const auto& c : output["coefficients"]) {
            csd_strings.push_back(c["csd"].get<std::string>());
        }

        auto max_len = static_cast<size_t>(0);
        for (const auto& s : csd_strings) {
            max_len = std::max(max_len, s.size());
        }
        auto max_power = static_cast<int>(max_len) - 1;

        std::vector<csd::MultiplierSpec> specs;
        for (size_t i = 0; i < csd_strings.size(); ++i) {
            auto raw = csd_strings[i];
            std::erase(raw, '.');
            while (raw.size() < max_len) {
                raw = "0" + raw;
            }
            specs.push_back({.name = "h" + std::to_string(i),
                             .csd = raw,
                             .input_width = input_width,
                             .max_power = max_power});
        }

        if (verilog_form == "transpose") {
            output["verilog"] = generate_transpose_form_verilog(specs, module_name);
        } else {
            auto raw = csd::generate_csd_multipliers(specs, module_name);
            output["verilog"] = fix_verilog_ports(raw);
        }
    }

    std::cout << output.dump(2) << '\n';
    return 0;
}
