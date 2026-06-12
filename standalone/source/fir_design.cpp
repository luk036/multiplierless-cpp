#include <algorithm>
#include <cmath>
#include <csd/csd_multiplier.hpp>
#include <ellalgo/arr.hpp>
#include <ellalgo/cutting_plane.hpp>
#include <ellalgo/ell.hpp>
#include <fstream>
#include <iostream>
#include <multiplierless/lowpass_oracle.hpp>
#include <multiplierless/lowpass_oracle_q.hpp>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

using json = nlohmann::json;

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

        std::vector<std::string> csd_strings;
        for (const auto& c : output["coefficients"]) {
            csd_strings.push_back(c["csd"].get<std::string>());
        }

        auto max_len = size_t(0);
        for (const auto& s : csd_strings) {
            max_len = std::max(max_len, s.size());
        }
        auto max_power = static_cast<int>(max_len) - 1;

        std::vector<csd::MultiplierSpec> specs;
        for (size_t i = 0; i < csd_strings.size(); ++i) {
            auto raw = csd_strings[i];
            raw.erase(std::remove(raw.begin(), raw.end(), '.'), raw.end());
            while (raw.size() < max_len) {
                raw = "0" + raw;
            }
            specs.push_back({"h" + std::to_string(i), raw, input_width, max_power});
        }

        output["verilog"] = csd::generate_csd_multipliers(specs, module_name);
    }

    std::cout << output.dump(2) << '\n';
    return 0;
}
