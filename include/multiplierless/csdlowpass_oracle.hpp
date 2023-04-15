#pragma once

#include <tuple>                       // for tuple
#include <type_traits>                 // for move
#include <xtensor/xlayout.hpp>         // for layout_type, layout_type::row...
#include <xtensor/xtensor_forward.hpp> // for xarray

#include "lowpass_oracle.hpp" // for LowpassOracle

class csdlowpass_oracle {
  using Arr = xt::xarray<double, xt::layout_type::row_major>;
  using Vec = std::valarray<double>;
  using ParallelCut = std::pair<Arr, Vec>;

private:
  Arr rcsd{};
  unsigned int _nnz;
  LowpassOracle _lowpass;

public:
  csdlowpass_oracle(unsigned int nnz, LowpassOracle &&lowpass)
      : _nnz(nnz), _lowpass(std::move(lowpass)) {}

  auto assess_q(const Arr &r, double &Spsq, bool retry)
      -> std::tuple<ParallelCut, bool, Arr, bool>;

  auto operator()(const Arr &r, double &Spsq, bool retry)
      -> std::tuple<ParallelCut, bool, Arr, bool> {
    return this->assess_q(r, Spsq, retry);
  }
};
