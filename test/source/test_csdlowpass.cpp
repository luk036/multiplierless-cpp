// -*- coding: utf-8 -*-
#include <doctest/doctest.h>

#include <cmath>
#include <complex>
#include <ellalgo/cutting_plane.hpp>
#include <ellalgo/ell.hpp>
#include <ellalgo/utility.hpp>
#include <limits>
#include <multiplierless/csdlowpass_oracle.hpp>
#include <tuple>
// #include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xview.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

template <int N> extern auto create_lowpass_case() -> std::tuple<lowpass_oracle, double>;

template <int N = 32> auto create_csdlowpass_case(int nnz = 8)
    -> std::tuple<csdlowpass_oracle, double> {
    auto [P, Spsq] = create_lowpass_case<N>();
    auto Pcsd = csdlowpass_oracle(nnz, std::move(P));
    return {Pcsd, Spsq};
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
    auto [P, t] = create_csdlowpass_case<N>(nnz);
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
