// -*- coding: utf-8 -*-
#include <doctest/doctest.h> // for ResultBuilder, CHECK

#include <ellalgo/cutting_plane.hpp>         // for cutting_plane_optim
#include <ellalgo/ell.hpp>                   // for Ell
#include <multiplierless/lowpass_oracle.hpp> // for filter_design_construct
#include <tuple>                             // for make_tuple, tuple
#include <type_traits>                       // for move, add_const<>::type
// #include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

// static filter_design_construct Fdc{};
auto create_lowpass_case(int N = 32) -> std::tuple<LowpassOracle, double> {
  auto Fdc = filter_design_construct(N);
  auto t = Fdc.Spsq;
  auto omega = LowpassOracle(std::move(Fdc));
  return {std::move(omega), t};
}

// ********************************************************************
// optimization
// ********************************************************************

auto run_lowpass(bool use_parallel_cut) {
  constexpr int N = 32;

  auto r0 = xt::zeros<double>({N}); // initial x0
  auto ellip = Ell<Arr>(40.0, r0);
  // auto omega = LowpassOracle(Fdc.Ap, Fdc.As, Fdc.Anr, Fdc.Lpsq, Fdc.Upsq);
  auto [omega, t] = create_lowpass_case(N);
  auto options = Options();

  options.max_iters = 50000;
  ellip.set_use_parallel_cut(use_parallel_cut);
  // options.tol = 1e-8;
  const auto [r, num_iters] = cutting_plane_optim(omega, ellip, t, options);
  // std::cout << "lowpass r: " << r << '\n';
  // auto Ustop = 20 * std::log10(std::sqrt(Spsq_new));
  // std::cout << "Min attenuation in the stopband is " << Ustop << " dB.\n";
  // CHECK(r[0] >= 0.0);
  return std::make_tuple(r.size() != 0U, num_iters);
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
