// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK

#include <ellalgo/cut_config.hpp>             // for Options
#include <ellalgo/cutting_plane.hpp>          // for cutting_plane_dc
#include <ellalgo/ell.hpp>                    // for ell
#include <multiplierless/lowpass_oracle.hpp>  // for filter_design_construct
#include <tuple>                              // for make_tuple, tuple
#include <type_traits>                        // for move, add_const<>::type
// #include <xtensor-blas/xlinalg.hpp>

// static filter_design_construct Fdc{};
auto create_lowpass_case(int N = 32) -> std::tuple<lowpass_oracle, double> {
    auto Fdc = filter_design_construct(N);
    auto t = Fdc.Spsq;
    auto P = lowpass_oracle(std::move(Fdc));
    return {std::move(P), t};
}

// ********************************************************************
// optimization
// ********************************************************************

auto run_lowpass(bool use_parallel_cut) {
    constexpr int N = 32;

    auto r0 = xt::zeros<double>({N});  // initial x0
    auto E = ell(40., r0);
    // auto P = lowpass_oracle(Fdc.Ap, Fdc.As, Fdc.Anr, Fdc.Lpsq, Fdc.Upsq);
    auto [P, t] = create_lowpass_case(N);
    auto options = Options();

    options.max_it = 50000;
    E.use_parallel_cut = use_parallel_cut;
    // options.tol = 1e-8;
    const auto [r, ell_info] = cutting_plane_dc(P, E, t, options);
    // std::cout << "lowpass r: " << r << '\n';
    // auto Ustop = 20 * std::log10(std::sqrt(Spsq_new));
    // std::cout << "Min attenuation in the stopband is " << Ustop << " dB.\n";
    // CHECK(r[0] >= 0.);
    return std::make_tuple(ell_info.feasible, ell_info.num_iters);
}

TEST_CASE("Lowpass Filter (w/ parallel cut)") {
    const auto [feasible, num_iters] = run_lowpass(true);
    CHECK(feasible);
    CHECK(num_iters <= 634);
}

// TEST_CASE("Lowpass Filter (w/o parallel cut)")
// {
//     const auto [feasible, num_iters] = run_lowpass(false);
//     CHECK(feasible);
//     CHECK(num_iters >= 7479);
// }
