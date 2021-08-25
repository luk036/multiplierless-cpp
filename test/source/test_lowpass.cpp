// -*- coding: utf-8 -*-
#include <doctest/doctest.h>

#include <cmath>
#include <complex>
#include <ellalgo/cutting_plane.hpp>
#include <ellalgo/ell.hpp>
#include <ellalgo/utility.hpp>
#include <limits>
#include <multiplierless/lowpass_oracle.hpp>
#include <tuple>
// #include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xview.hpp>

#ifndef M_PI
#    define M_PI 3.14159265358979323846264338327950288
#endif

using Arr = xt::xarray<double, xt::layout_type::row_major>;
// static const auto M_PI = std::acos(-1);

// Modified from CVX code by Almir Mutapcic in 2006.
// Adapted in 2010 for impulse response peak-minimization by convex iteration by
// Christine Law.
//
// "FIR Filter Design via Spectral Factorization and Convex Optimization"
// by S.-P. Wu, S. Boyd, and L. Vandenberghe
//
// Designs an FIR lowpass filter using spectral factorization method with
// constraint on maximum passband ripple and stopband attenuation:
//
//   minimize   max |H(w)|                      for w in stopband
//       s.t.   1/delta <= |H(w)| <= delta      for w in passband
//
// We change variables via spectral factorization method and get:
//
//   minimize   max R(w)                          for w in stopband
//       s.t.   (1/delta)**2 <= R(w) <= delta**2  for w in passband
//              R(w) >= 0                         for all w
//
// where R(w) is squared magnitude frequency response
// (and Fourier transform of autocorrelation coefficients r).
// Variables are coeffients r and G = hh' where h is impulse response.
// delta is allowed passband ripple.
// This is a convex problem (can be formulated as an SDP after sampling).

// rand('twister',sum(100*clock))
// randn('state',sum(100*clock))

// *********************************************************************
// filter specs (for a low-pass filter)
// *********************************************************************
// number of FIR coefficients (including zeroth)
struct filter_design_construct {
    int N;
    Arr Ap;
    Arr As;
    Arr Anr;
    double Lpsq;
    double Upsq;
    double Spsq;

    filter_design_construct(int argN = 32) : N(argN) {
        const auto wpass = 0.12 * M_PI;  // end of passband
        const auto wstop = 0.20 * M_PI;  // start of stopband
        const auto delta0_wpass = 0.125;
        const auto delta0_wstop = 0.125;
        // maximum passband ripple in dB (+/- around 0 dB)
        const auto delta = 20 * std::log10(1 + delta0_wpass);
        // stopband attenuation desired in dB
        const auto delta2 = 20 * std::log10(delta0_wstop);

        // *********************************************************************
        // optimization parameters
        // *********************************************************************
        // rule-of-thumb discretization (from Cheney's Approximation Theory)
        const auto m = 15 * this->N;
        const auto w = Arr{xt::linspace<double>(0, M_PI, m)};  // omega

        // passband 0 <= w <= w_pass
        const auto Lp = std::pow(10, -delta / 20);
        const auto Up = std::pow(10, +delta / 20);

        // A is the matrix used to compute the power spectrum
        // A(w,:) = [1 2*cos(w) 2*cos(2*w) ... 2*cos((this->N-1)*w)]

        // Arr An = 2 * xt::cos(xt::linalg::outer(w, xt::arange(1, this->N)));
        auto An = Arr(xt::zeros<double>({m, this->N - 1}));
        for (auto i = 0; i != m; ++i) {
            for (auto j = 0; j != this->N - 1; ++j) {
                An(i, j) = 2. * std::cos(w(i) * (j + 1));
            }
        }
        Arr A = xt::concatenate(xt::xtuple(xt::ones<double>({m, 1}), An), 1);

        const auto ind_p = xt::where(w <= wpass)[0];  // passband
        this->Ap = Arr{xt::view(A, xt::range(0, ind_p.size()), xt::all())};

        // stopband (w_stop <= w)
        auto ind_s = xt::where(wstop <= w)[0];  // stopband
        const auto Sp = std::pow(10, delta2 / 20);

        using xt::placeholders::_;
        this->As = xt::view(A, xt::range(ind_s[0], _), xt::all());

        // remove redundant contraints
        // ind_nr = setdiff(1:m,ind_p)   // fullband less passband
        // ind_nr = setdiff(ind_nr, ind_s) // luk: for making parallel cut
        // auto ind_nr = np.setdiff1d(xt::arange(m), ind_p);
        // auto ind_nr = np.setdiff1d(ind_nr, ind_s);
        auto ind_beg = ind_p[ind_p.size() - 1];
        auto ind_end = ind_s[0];
        this->Anr = xt::view(A, xt::range(ind_beg + 1, ind_end), xt::all());

        this->Lpsq = Lp * Lp;
        this->Upsq = Up * Up;
        this->Spsq = Sp * Sp;
    }
};

// static filter_design_construct Fdc{};
template <int N = 32> auto create_lowpass_case() -> std::tuple<lowpass_oracle, double> {
    static auto Fdc = filter_design_construct(N);
    auto P = lowpass_oracle(Fdc.Ap, Fdc.As, Fdc.Anr, Fdc.Lpsq, Fdc.Upsq);
    return {std::move(P), Fdc.Spsq};
}

// ********************************************************************
// optimization
// ********************************************************************

auto run_lowpass(bool use_parallel_cut) {
    constexpr int N = 32;

    auto r0 = zeros({N});  // initial x0
    auto E = ell(40., r0);
    // auto P = lowpass_oracle(Fdc.Ap, Fdc.As, Fdc.Anr, Fdc.Lpsq, Fdc.Upsq);
    auto [P, t] = create_lowpass_case<N>();
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
