#include <iosfwd> // for string
#include <multiplierless/csdlowpass_oracle.hpp>
#include <tuple>                  // for tuple_element<>::type
#include <xtensor/xarray.hpp>     // for xarray_container
#include <xtensor/xcontainer.hpp> // for xcontainer
#include <xtensor/xfunction.hpp>  // for xfunction
#include <xtensor/xmath.hpp>      // for sum
#include <xtensor/xoperation.hpp> // for xfunction_type_t, opera...
#include <xtensor/xreducer.hpp>   // for xreducer, xreducer<>::c...
#include <xtensor/xsemantic.hpp>  // for xsemantic_base

#include "multiplierless/lowpass_oracle.hpp" // for LowpassOracle

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using Vec = std::valarray<double>;
using ParallelCut = std::pair<Arr, Vec>;

extern auto to_csdfixed(double num, unsigned int) -> std::string;
extern auto to_decimal(const std::string &csd_str) -> double;

// extern auto to_csdfixed();
// extern auto to_decimal();
extern auto inverse_spectral_fact(const Arr &r) -> Arr;
extern auto spectral_fact(const Arr &r) -> Arr;

auto csdlowpass_oracle::assess_q(const Arr &r, double &Spsq, bool retry)
    -> std::tuple<ParallelCut, bool, Arr, bool> {
  if (!retry) { // retry due to no effect in the previous cut
    this->_lowpass.retry = false;
    auto [cut, shrunk] = this->_lowpass(r, Spsq);
    if (!shrunk) {
      return {cut, shrunk, r, true};
    }
    auto h = spectral_fact(r);
    auto hcsd = Arr(h.shape());
    for (auto i = 0U; i != h.size(); ++i) {
      hcsd(i) = to_decimal(to_csdfixed(h(i), this->_nnz));
    }
    this->rcsd = inverse_spectral_fact(hcsd);
  }
  auto [cut, shrunk] = this->_lowpass(this->rcsd, Spsq);
  auto &[gc, hc] = cut;
  // no more alternative cuts?
  hc += xt::sum(gc * (this->rcsd - r))();
  auto more_alt = this->_lowpass.more_alt && !retry;
  return {cut, shrunk, this->rcsd, more_alt};
}
