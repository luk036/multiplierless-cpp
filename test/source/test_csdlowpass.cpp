// -*- coding: utf-8 -*-
#include <ellalgo/cutting_plane.hpp>            // for cutting_plane_q
#include <ellalgo/ell.hpp>                      // for Ell
#include <multiplierless/csdlowpass_oracle.hpp> // for LowpassOracleQ
#include <tuple>                                // for make_tuple, tuple
#include <type_traits>                          // for move, add_const<>::type
#include <xtensor/xlayout.hpp>                  // for layout_type, layout_...
#include <xtensor/xtensor_forward.hpp>          // for xarray

class LowpassOracle;

using Arr = xt::xarray<double, xt::layout_type::row_major>;

extern auto create_lowpass_case(int N) -> std::tuple<LowpassOracle, double>;

auto create_csdlowpass_case(int N = 32, int nnz = 8)
    -> std::tuple<LowpassOracleQ, double> {
  auto [omega, Spsq] = create_lowpass_case(N);
  auto Pcsd = LowpassOracleQ(nnz, std::move(omega));
  return {std::move(Pcsd), Spsq};
}

// ********************************************************************
// optimization
// ********************************************************************

auto run_csdlowpass(bool use_parallel_cut) {
  constexpr int N = 32;
  const int nnz = 7;

  auto r0 = xt::zeros<double>({N}); // initial x0
  auto ellip = Ell<Arr>(40.0, r0);
  // auto omega = LowpassOracleQ(Fdc.Ap, Fdc.As, Fdc.Anr, Fdc.Lpsq, Fdc.Upsq);
  auto [omega, t] = create_csdlowpass_case(N, nnz);
  auto options = Options();

  options.max_iter = 50000;
  ellip.set_use_parallel_cut(use_parallel_cut);
  // options.tol = 1e-8;

  // auto t = Fdc.Spsq;
  const auto [r, num_iters] = cutting_plane_q(omega, ellip, t, options);
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
