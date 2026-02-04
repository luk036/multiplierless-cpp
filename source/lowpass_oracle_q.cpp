#include <iosfwd>  // for string
#include <multiplierless/lowpass_oracle_q.hpp>
#include <tuple>                   // for tuple_element<>::type
#include <xtensor/xarray.hpp>      // for xarray_container
#include <xtensor/xcontainer.hpp>  // for xcontainer
#include <xtensor/xfunction.hpp>   // for xfunction
#include <xtensor/xmath.hpp>       // for sum
#include <xtensor/xoperation.hpp>  // for xfunction_type_t, opera...
#include <xtensor/xreducer.hpp>    // for xreducer, xreducer<>::c...
#include <xtensor/xsemantic.hpp>   // for xsemantic_base

#include "multiplierless/lowpass_oracle.hpp"  // for LowpassOracle

using Arr = xt::xarray<double>;
using Vec = std::valarray<double>;
using ParallelCut = std::pair<Arr, Vec>;

extern auto to_csdnnz(double num, unsigned int) -> std::string;
extern auto to_decimal(const std::string &csd_str) -> double;

extern auto inverse_spectral_fact(const Arr &r) -> Arr;
extern auto spectral_fact(const Arr &r) -> Arr;

/**
 * The function assess_optim_q assesses the optimal value of q for a lowpass filter based on the
 * input signal and previous cuts.
 *
 * @param[in] r The parameter `r` is the top-half of the auto-correlation coefficients. It
 * represents the desired auto-correlation of the impulse response. It should be passed in as a
 * column vector.
 * @param[in] Spsq Spsq is a reference to a double variable. It is used to store the sum of squared
 * values of the input signal.
 * @param[in] retry A boolean flag indicating whether the function is being retried due to no effect
 * in the previous cut.
 *
 * @return The function `assess_optim_q` returns a tuple containing the following elements:
 */
auto LowpassOracleQ::assess_optim_q(const Arr &r, double &Spsq, bool retry)
    -> std::tuple<ParallelCut, bool, Arr, bool> {
    if (!retry) {  // retry due to no effect in the previous cut
        // this->_lowpass.retry = false;
        auto [cut, shrunk] = this->_lowpass(r, Spsq);
        if (!shrunk) {
            return {cut, shrunk, r, true};
        }
        auto h = spectral_fact(r);
        auto hcsd = Arr(h.shape());
        for (auto i = 0U; i != h.size(); ++i) {
            hcsd(i) = to_decimal(to_csdnnz(h(i), this->_nnz));
        }
        this->rcsd = inverse_spectral_fact(hcsd);
        this->_num_retries = 0;  // reset to zero
    } else {
        this->_num_retries += 1;
    }
    auto [cut, shrunk] = this->_lowpass(this->rcsd, Spsq);
    auto &[gc, hc] = cut;
    hc += xt::sum(gc * (this->rcsd - r))();
    // auto more_alt = this->_lowpass.more_alt && !retry;
    return {cut, shrunk, this->rcsd, this->_num_retries < 15};
}
