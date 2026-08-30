// -*- coding: utf-8 -*-
#include <ellalgo/arr.hpp>
#include <ellalgo/cutting_plane.hpp>            // for cutting_plane_optim_q
#include <ellalgo/ell.hpp>                      // for Ell
#include <multiplierless/lowpass_oracle_q.hpp>  // for LowpassOracleQ, create_csdlowpass_case
#include <tuple>                                // for make_tuple, tuple

// ********************************************************************
// optimization
// ********************************************************************

/**
 * @brief Run the CSD-quantised lowpass optimisation.
 * @param[in] use_parallel_cut Whether to enable parallel cutting-plane.
 * @return (feasible, iteration_count) tuple.
 */
auto run_csdlowpass(bool use_parallel_cut) {
    constexpr int N = 32;
    const int nnz = 7;

    auto r0 = zeros(N);  // initial x0
    auto ellip = Ell<Arr>(40.0, r0);
    // auto omega = LowpassOracleQ(Fdc.Ap, Fdc.As, Fdc.Anr, Fdc.Lpsq, Fdc.Upsq);
    auto [omega, t] = create_csdlowpass_case(N, nnz);
    auto options = Options();

    options.max_iters = 50000;
    options.tolerance = 1e-14;
    ellip.set_use_parallel_cut(use_parallel_cut);

    // auto t = Fdc.Spsq;
    const auto [r, num_iters] = cutting_plane_optim_q(omega, ellip, t, options);
    // std::cout << "lowpass r: " << r << '\n';
    // auto Ustop = 20 * std::log10(std::sqrt(Spsq_new));
    // std::cout << "Min attenuation in the stopband is " << Ustop << " dB.\n";
    // CHECK(r[0] >= 0.0);
    return std::make_tuple(r.size() != 0U, num_iters);
}

// TEST_CASE("CSD Lowpass Filter (w/ parallel cut)") {
//     const auto [feasible, num_iters] = run_csdlowpass(true);
//     CHECK(feasible);
//     CHECK(num_iters <= 1136);
// }

// TEST_CASE("Lowpass Filter (w/o parallel cut)")
// {
//     const auto [feasible, num_iters] = run_lowpass(false);
//     CHECK(feasible);
//     CHECK(num_iters >= 7479);
// }
