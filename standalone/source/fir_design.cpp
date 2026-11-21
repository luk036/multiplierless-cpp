#include <algorithm>
#include <cmath>
#include <csd/csd_multiplier.hpp>
#include <ellalgo/arr.hpp>
#include <ellalgo/cutting_plane.hpp>
#include <ellalgo/ell.hpp>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <multiplierless/csd_quantize.hpp>
#include <multiplierless/lowpass_oracle.hpp>
#include <multiplierless/lowpass_oracle_q.hpp>
#include <multiplierless/spectral_fact.hpp>
#include <multiplierless/verilog_utils.hpp>
#include <nlohmann/json.hpp>
#include <regex>
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
        for (size_t ci = 0; ci < csd_list.size(); ++ci) {
            auto const& csd = csd_list[ci];
            auto const n = csd.size();
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 2; j <= n; ++j) {
                    auto sub = csd.substr(i, j - i);
                    if (count_nnz(sub) >= 2) {
                        patterns[sub].emplace_back(static_cast<int>(ci), static_cast<int>(i));
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
        // Summing N products can grow by log2(N) bits beyond a single product.
        auto const headroom
            = (N > 1) ? static_cast<int>(std::ceil(std::log2(static_cast<double>(N)))) : 0;
        auto const output_width = input_width + max_power + headroom;

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
        std::set<int, std::greater<>> all_powers;
        for (auto const& spec : coeffs) {
            for (size_t i = 0; i < spec.csd.size(); ++i) {
                if (spec.csd[i] != '0') {
                    all_powers.insert(max_power - static_cast<int>(i));
                }
            }
        }

        // Cross-CSE pattern detection
        std::vector<std::string> csd_strings;
        csd_strings.reserve(static_cast<size_t>(N));
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
                cse_base_pos = std::min(occ_pair.second, cse_base_pos);
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
            auto const& spec = coeffs[static_cast<size_t>(coeff_idx)];
            bool has_cse = !best_pattern.empty() && (cse_coeffs.contains(coeff_idx));
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

    // Strategy registry: verilog form -> generator (transpose / direct)
    using VerilogGenerator
        = std::function<std::string(const std::vector<csd::MultiplierSpec>&, const std::string&)>;
    const std::map<std::string, VerilogGenerator> verilog_generators{
        {"transpose",
         [](const auto& specs, const std::string& module_name) {
             return generate_transpose_form_verilog(specs, module_name);
         }},
        {"direct",
         [](const auto& specs, const std::string& module_name) {
             std::vector<std::string> names;
             names.reserve(specs.size());
             for (auto const& spec : specs) names.push_back(spec.name);
             return fix_direct_verilog(csd::generate_csd_multipliers(specs, module_name), names);
         }},
    };
}  // anonymous namespace

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
namespace {

    struct CliArgs {
        std::string spec_path;
        std::string output_path;
        std::string verilog_path;
        bool show_help = false;
    };

    auto usage(std::ostream& os, const char* prog) -> void {
        os << "Usage: " << prog << " <filter_spec.json> [options]\n"
           << "  -o, --output FILE       write result JSON to FILE (default: stdout)\n"
           << "      --verilog-out FILE  write generated Verilog to FILE\n"
           << "  -h, --help              show this help\n";
    }

    auto parse_args(int argc, char** argv) -> CliArgs {
        CliArgs args;
        for (int i = 1; i < argc; ++i) {
            std::string const arg = argv[i];
            if (arg == "-h" || arg == "--help") {
                args.show_help = true;
            } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
                args.output_path = argv[++i];
            } else if (arg == "--verilog-out" && i + 1 < argc) {
                args.verilog_path = argv[++i];
            } else if (!arg.empty() && arg[0] == '-') {
                throw std::invalid_argument("unknown option: " + arg);
            } else if (args.spec_path.empty()) {
                args.spec_path = arg;
            } else {
                throw std::invalid_argument("unexpected argument: " + arg);
            }
        }
        return args;
    }

    auto write_text(const std::string& path, const std::string& text) -> void {
        std::ofstream out(path, std::ios::binary);
        if (!out) {
            throw std::runtime_error("cannot write " + path);
        }
        out << text;
        if (!out) {
            throw std::runtime_error("failed while writing " + path);
        }
    }

}  // namespace

int main(int argc, char** argv) {
    try {
        auto const args = parse_args(argc, argv);
        if (args.show_help) {
            usage(std::cout, argv[0]);
            return 0;
        }
        if (args.spec_path.empty()) {
            usage(std::cerr, argv[0]);
            return 1;
        }

        std::ifstream ifs(args.spec_path);
        if (!ifs) {
            std::cerr << "error: cannot open " << args.spec_path << '\n';
            return 1;
        }

        json spec;
        try {
            ifs >> spec;
        } catch (const json::parse_error& e) {
            std::cerr << "error: invalid JSON in " << args.spec_path << ": " << e.what() << '\n';
            return 1;
        }
        if (!spec.is_object()) {
            std::cerr << "error: filter specification must be a JSON object\n";
            return 1;
        }

        auto filter_order = spec.value("filter_order", 32);
        auto csd_nnz = spec.value("csd_nnz", 7U);
        if (filter_order < 1) {
            throw std::invalid_argument("filter_order must be a positive integer");
        }
        if (csd_nnz < 1 || csd_nnz > 16) {
            throw std::invalid_argument("csd_nnz must be between 1 and 16");
        }

        auto spectral_method_str = spec.value("spectral_method", std::string("fft"));
        auto method = spectral_method::fft;
        if (spectral_method_str == "fft") {
            method = spectral_method::fft;
        } else if (spectral_method_str == "root") {
            method = spectral_method::root;
        } else {
            throw std::invalid_argument("unknown spectral_method: " + spectral_method_str);
        }
        auto root_tol = spec.value("root_tolerance", 1e-8);
        if (root_tol <= 0.0) {
            throw std::invalid_argument("root_tolerance must be positive");
        }

        auto wpass_norm = spec.value("passband_edge", 0.12);
        auto wstop_norm = spec.value("stopband_edge", 0.20);
        auto passband_ripple = spec.value("passband_ripple", 0.125);
        auto stopband_attn = spec.value("stopband_attenuation", 0.125);
        auto disc_factor = spec.value("discretization_factor", 15);
        auto max_iters = spec.value("max_iters", static_cast<size_t>(50000));
        auto tolerance = spec.value("tolerance", 1e-14);
        auto ell_radius = spec.value("ellipsoid_radius", 40.0);
        auto parallel_cut = spec.value("parallel_cut", true);
        auto gen_verilog = spec.contains("verilog");
        if (gen_verilog && !spec["verilog"].is_object()) {
            throw std::invalid_argument("'verilog' must be a JSON object");
        }

        auto Fdc = filter_design_construct(filter_order, wpass_norm, wstop_norm, passband_ripple,
                                           stopband_attn, disc_factor);
        auto Spsq = Fdc.Spsq;
        auto omega = LowpassOracleQ(csd_nnz, LowpassOracle(std::move(Fdc)));

        auto r0 = zeros(static_cast<size_t>(filter_order));
        auto ellip = Ell<Arr>(ell_radius, r0);
        ellip.set_use_parallel_cut(parallel_cut);

        Options options(max_iters, tolerance);
        auto [r, num_iters] = cutting_plane_optim_q(omega, ellip, Spsq, options);
        if (r.size() == 0) {
            std::cerr << "error: no feasible solution after " << num_iters << " iterations\n";
            return 1;
        }

        Arr h = spectral_fact(r, method, root_tol);

        json output;
        output["filter_order"] = filter_order;
        output["csd_nnz"] = csd_nnz;
        output["iterations"] = num_iters;
        output["spectral_method"] = spectral_method_str;
        output["coefficients"] = json::array();
        for (size_t i = 0; i < h.size(); ++i) {
            auto const parts = csd_quantize_parts(h(i), csd_nnz);
            output["coefficients"].push_back(
                {{"index", i}, {"value", parts.value}, {"csd", format_csd(parts.digits)}});
        }

        std::string verilog_text;
        if (gen_verilog) {
            auto const& vl = spec["verilog"];
            auto input_width = vl.value("input_width", 16);
            auto module_name = vl.value("module_name", std::string("fir_filter"));
            auto verilog_form = vl.value("form", std::string("transpose"));
            if (input_width < 1) {
                throw std::invalid_argument("verilog.input_width must be a positive integer");
            }
            if (!std::regex_match(module_name, std::regex("[A-Za-z_][A-Za-z0-9_]*"))) {
                throw std::invalid_argument("invalid verilog.module_name: " + module_name);
            }
            auto const generator = verilog_generators.find(verilog_form);
            if (generator == verilog_generators.end()) {
                throw std::invalid_argument("unknown verilog form: " + verilog_form);
            }

            std::vector<std::string> csd_strings;
            csd_strings.reserve(output["coefficients"].size());
            for (auto const& c : output["coefficients"]) {
                csd_strings.push_back(c["csd"].get<std::string>());
            }
            auto const fixed = csd_coeffs_to_fixed(csd_strings);
            std::vector<csd::MultiplierSpec> specs;
            specs.reserve(fixed.digits.size());
            for (size_t i = 0; i < fixed.digits.size(); ++i) {
                specs.push_back({.name = "h" + std::to_string(i),
                                 .csd = fixed.digits[i],
                                 .input_width = input_width,
                                 .max_power = fixed.max_power});
            }
            verilog_text = generator->second(specs, module_name);
            output["verilog"] = verilog_text;
        }

        if (!args.verilog_path.empty()) {
            if (!gen_verilog) {
                throw std::invalid_argument("--verilog-out requires a 'verilog' section");
            }
            write_text(args.verilog_path, verilog_text);
        }

        auto const text = output.dump(2) + '\n';
        if (!args.output_path.empty()) {
            write_text(args.output_path, text);
        } else {
            std::cout << text;
        }
        return 0;
    } catch (const json::exception& e) {
        std::cerr << "error: invalid filter specification: " << e.what() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
    }
    return 1;
}
