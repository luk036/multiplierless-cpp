#pragma once

#include <tuple>                        // for tuple
#include <type_traits>                  // for move
#include <xtensor/xlayout.hpp>          // for layout_type, layout_type::row...
#include <xtensor/xtensor_forward.hpp>  // for xarray

#include "lowpass_oracle.hpp"  // for LowpassOracle

/* The `LowpassOracleQ` class is a wrapper class that provides an interface for optimizing a lowpass filter. It takes an instance of the `LowpassOracle` class as a parameter and stores it as a member variable `_lowpass`. */
class LowpassOracleQ {
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Vec = std::valarray<double>;
    using ParallelCut = std::pair<Arr, Vec>;

    Arr rcsd{};
    unsigned int _nnz;
    LowpassOracle _lowpass;

  public:
    /* The `LowpassOracleQ` constructor is initializing a `LowpassOracleQ` object. It takes two parameters: `nnz` of type `unsigned int` and `lowpass` of type `LowpassOracle&&` (an rvalue reference to a `LowpassOracle` object). */
    LowpassOracleQ(unsigned int nnz, LowpassOracle &&lowpass)
        : _nnz(nnz), _lowpass(std::move(lowpass)) {}

    /* The `assess_optim_q` function is a member function of the `LowpassOracleQ` class. It takes three parameters: `r` of type `const Arr&` (a reference to a constant `Arr` object), `Spsq` of type `double&` (a reference to a `double`), and `retry` of type `bool`. */
    auto assess_optim_q(const Arr &r, double &Spsq, bool retry)
        -> std::tuple<ParallelCut, bool, Arr, bool>;

    auto operator()(const Arr &r, double &Spsq, bool retry)
        -> std::tuple<ParallelCut, bool, Arr, bool> {
        return this->assess_optim_q(r, Spsq, retry);
    }
};
