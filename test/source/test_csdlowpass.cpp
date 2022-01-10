// -*- coding: utf-8 -*-
#include <ellalgo/cut_config.hpp>                // for Options
#include <ellalgo/cutting_plane.hpp>             // for cutting_plane_q
#include <ellalgo/ell.hpp>                       // for ell
#include <ellalgo/utility.hpp>                   // for zeros
#include <multiplierless/csdlowpass_oracle.hpp>  // for csdlowpass_oracle
#include <tuple>                                 // for make_tuple, tuple
#include <type_traits>                           // for move, add_const<>::type
#include <xtensor/xlayout.hpp>                   // for layout_type, layout_...
#include <xtensor/xtensor_forward.hpp>           // for xarray

class lowpass_oracle;

using Arr = xt::xarray<double, xt::layout_type::row_major>;

extern auto create_lowpass_case(int N) -> std::tuple<lowpass_oracle, double>;

auto create_csdlowpass_case(int N = 32, int nnz = 8) -> std::tuple<csdlowpass_oracle, double> {
    auto [P, Spsq] = create_lowpass_case(N);
    auto Pcsd = csdlowpass_oracle(nnz, std::move(P));
    return {std::move(Pcsd), Spsq};
}

// ********************************************************************
// optimization
// ********************************************************************

auto run_csdlowpass(bool use_parallel_cut) {
    constexpr int N = 32;
    const int nnz = 7;

    auto r0 = zeros({N});  // initial x0
    auto E = ell(40., r0);
    // auto P = csdlowpass_oracle(Fdc.Ap, Fdc.As, Fdc.Anr, Fdc.Lpsq, Fdc.Upsq);
    auto [P, t] = create_csdlowpass_case(N, nnz);
    auto options = Options();

    options.max_it = 50000;
    E.use_parallel_cut = use_parallel_cut;
    // options.tol = 1e-8;

    // auto t = Fdc.Spsq;
    const auto [r, ell_info] = cutting_plane_q(P, E, t, options);
    // std::cout << "lowpass r: " << r << '\n';
    // auto Ustop = 20 * std::log10(std::sqrt(Spsq_new));
    // std::cout << "Min attenuation in the stopband is " << Ustop << " dB.\n";
    // CHECK(r[0] >= 0.);
    return std::make_tuple(ell_info.feasible, ell_info.num_iters);
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
