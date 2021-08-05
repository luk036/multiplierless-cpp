#pragma once

#include "lowpass_oracle.hpp"

class csdlowpass_oracle {
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using ParallelCut = std::tuple<Arr, Arr>;

  private:
    Arr rcsd{};
    unsigned int _nnz;
    lowpass_oracle _lowpass;

  public:
    csdlowpass_oracle(unsigned int nnz, lowpass_oracle&& lowpass) : 
        _nnz(nnz), _lowpass(std::move(lowpass)) {}

    auto operator()(const Arr& r, double& Spsq, bool retry)
        -> std::tuple<ParallelCut, bool, Arr, bool>;
};