#include <multiplierless/csdlowpass_oracle.hpp>
#include <xtensor/xview.hpp>
#include <string>
#include <string_view>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using ParallelCut = std::tuple<Arr, Arr>;

extern auto to_csdfixed(double num, unsigned int ) -> std::string;
extern auto to_decimal(std::string_view csd_str) -> double;

extern auto to_csdfixed();
extern auto to_decimal();
extern auto inverse_spectral_fact(const Arr& r) -> Arr;
extern auto spectral_fact(const Arr& r) -> Arr;

auto csdlowpass_oracle::operator()(const Arr& r, double& Spsq, bool retry)
    -> std::tuple<ParallelCut, bool, Arr, bool>
{
    if (!retry) {  // retry due to no effect in the previous cut
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
    auto& [gc, hc] = cut;
    // no more alternative cuts?
    hc += xt::sum(gc * (this->rcsd - r))();
    auto more_alt = this->_lowpass.more_alt && !retry;
    return {cut, shrunk, this->rcsd, more_alt};
}