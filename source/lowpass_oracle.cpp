#include <stddef.h> // for size_t

#include <cmath>                             // for pow, log10, M_PI, cos
#include <multiplierless/lowpass_oracle.hpp> // for LowpassOracle, filter_...
#include <tuple>                             // for tuple
#include <type_traits>                       // for move
#include <vector>                            // for vector, vector<>::size_...
#include <xtensor/xaccessible.hpp>           // for xconst_accessible, xacc...
#include <xtensor/xarray.hpp>                // for xarray_container
#include <xtensor/xbroadcast.hpp>            // for xbroadcast
#include <xtensor/xbuilder.hpp>              // for zeros, concatenate, lin...
#include <xtensor/xcontainer.hpp>            // for xcontainer<>::inner_sha...
#include <xtensor/xexception.hpp>            // for throw_concatenate_error
#include <xtensor/xgenerator.hpp>            // for xgenerator
#include <xtensor/xlayout.hpp>               // for layout_type, layout_typ...
#include <xtensor/xmath.hpp>                 // for sum
#include <xtensor/xoperation.hpp>            // for xfunction_type_t, opera...
#include <xtensor/xreducer.hpp>              // for xreducer
#include <xtensor/xslice.hpp>                // for all, range, xtuph, _
#include <xtensor/xtensor_forward.hpp>       // for xarray
#include <xtensor/xutils.hpp>                // for accumulate
#include <xtensor/xview.hpp>                 // for xview, view
                                             //
// #include <limits>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using Vec = std::valarray<double>;
using ParallelCut = std::pair<Arr, Vec>;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

filter_design_construct::filter_design_construct(int argN) : N(argN) {
    const auto wpass = 0.12 * M_PI; // end of passband
    const auto wstop = 0.20 * M_PI; // start of stopband
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
    const auto w = Arr{xt::linspace<double>(0, M_PI, size_t(m))}; // omega
    // passband 0 <= w <= w_pass
    const auto Lp = std::pow(10, -delta / 20);
    const auto Up = std::pow(10, +delta / 20);
    // A is the matrix used to compute the power spectrum
    // A(w,:) = [1 2*cos(w) 2*cos(2*w) ... 2*cos((this->N-1)*w)]
    // Arr An = 2 * xt::cos(xt::linalg::outer(w, xt::arange(1, this->N)));
    Arr An = xt::zeros<double>({m, this->N - 1});
    for (auto i = 0; i != m; ++i) {
        for (auto j = 0; j != this->N - 1; ++j) {
            An(i, j) = 2.0 * std::cos(w(i) * (j + 1));
        }
    }
    Arr A = xt::concatenate(xt::xtuple(xt::ones<double>({m, 1}), An), 1);
    const auto ind_p = xt::where(w <= wpass)[0]; // passband
    this->Ap = xt::view(A, xt::range(0, ind_p.size()), xt::all());
    // stopband (w_stop <= w)
    auto ind_s = xt::where(wstop <= w)[0]; // stopband
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

/*!
 * @brief
 *
 * @param[in] x
 * @param[in] Spsq
 * @return auto
 */
auto LowpassOracle::assess_optim(const Arr &x, double &Spsq)
    -> std::tuple<ParallelCut, bool> {
    this->more_alt = true;

    // 1.0 nonnegative-real constraint
    // case 1,
    if (x[0] < 0) {
        Arr g = xt::zeros<double>(x.shape());
        g[0] = -1.;
        auto f = Vec{-x[0]};
        return {{std::move(g), std::move(f)}, false};
    }

    // case 2,
    // 2.0 passband constraints
    auto N = this->_Fdc.Ap.shape()[0];
    // for (k in chain(range(i_As, N), range(i_As))) {

    this->retry = false; // ???

    auto k = this->_i_Ap;
    for (auto i = 0U; i != N; ++i, ++k) {
        if (k == N) {
            k = 0; // round robin
        }
        auto v = xt::sum(xt::view(this->_Fdc.Ap, k, xt::all()) * x)();
        if (v > this->_Fdc.Upsq) {
            // f = v - Upsq;
            Arr g = xt::view(this->_Fdc.Ap, k, xt::all());
            Vec f{v - this->_Fdc.Upsq, v - this->_Fdc.Lpsq};
            this->_i_Ap = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
        if (v < this->_Fdc.Lpsq) {
            // f = Lpsq - v;
            Arr g = -xt::view(this->_Fdc.Ap, k, xt::all());
            Vec f{-v + this->_Fdc.Lpsq, -v + this->_Fdc.Upsq};
            this->_i_Ap = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
    }

    // case 3,
    // 3.0 stopband constraint
    N = this->_Fdc.As.shape()[0];
    // Arr w = xt::zeros<double>({N});
    auto fmax = -1.e100; // std::numeric_limits<double>::min()
    size_t imax = 0U;
    // for (k in chain(range(i_As, N), range(i_As))) {
    k = this->_i_As;
    for (auto i = 0U; i != N; ++i, ++k) {
        if (k == N) {
            k = 0; // round robin
        }
        auto v = xt::sum(xt::view(this->_Fdc.As, k, xt::all()) * x)();
        if (v > Spsq) {
            // f = v - Spsq;
            Arr g = xt::view(this->_Fdc.As, k, xt::all());
            // f = (v - Spsq, v);
            Vec f{v - Spsq, v};
            this->_i_As = k + 1; // k or k+1
            return {{std::move(g), std::move(f)}, false};
        }
        if (v < 0) {
            // f = v - Spsq;
            Arr g = -xt::view(this->_Fdc.As, k, xt::all());
            Vec f{-v, -v + Spsq};
            this->_i_As = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
        if (v > fmax) {
            fmax = v;
            imax = k;
        }
    }

    // case 4,
    // 1.0 nonnegative-real constraint
    N = this->_Fdc.Anr.shape()[0];
    k = this->_i_Anr;
    for (auto i = 0U; i != N; ++i, ++k) {
        if (k == N) {
            k = 0; // round robin
        }
        auto v = xt::sum(xt::view(this->_Fdc.Anr, k, xt::all()) * x)();
        if (v < 0.0) {
            Vec f{-v};
            Arr g = -xt::view(this->_Fdc.Anr, k, xt::all());
            this->_i_Anr = k + 1;
            return {{std::move(g), std::move(f)}, false};
        }
    }

    this->more_alt = false;

    // Begin objective function
    // Spsq, imax = w.max(), w.argmax(); // update best so far Spsq
    Spsq = fmax;
    Vec f{0.0, fmax}; // ???
    // f = 0
    Arr g = xt::view(this->_Fdc.As, imax, xt::all());
    return {{std::move(g), std::move(f)}, true};
}
