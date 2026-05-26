#include <ellalgo/arr.hpp>
#include <multiplierless/lowpass_oracle_q.hpp>

using Vec = std::valarray<double>;
using ParallelCut = std::pair<Arr, Vec>;

extern auto to_csdnnz_fast(double num, unsigned int) -> double;

extern auto inverse_spectral_fact(const Arr& r) -> Arr;
extern auto spectral_fact(const Arr& r) -> Arr;

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
auto LowpassOracleQ::assess_optim_q(const Arr& r, double& Spsq, bool retry)
    -> std::tuple<ParallelCut, bool, Arr, bool> {
    if (!retry) {  // retry due to no effect in the previous cut
        auto [cut, shrunk] = this->_lowpass(r, Spsq);
        if (!shrunk) {
            return {cut, shrunk, r, true};
        }
        auto h = spectral_fact(r);
        auto hcsd = Arr(h.size());
        for (auto i = 0U; i != h.size(); ++i) {
            hcsd(i) = to_csdnnz_fast(h(i), this->_nnz);
        }
        this->rcsd = inverse_spectral_fact(hcsd);
        this->_num_retries = 0;  // reset to zero
    } else {
        this->_num_retries += 1;
    }
    auto [cut, shrunk] = this->_lowpass(this->rcsd, Spsq);
    auto& [gc, hc] = cut;
    auto dot_val = 0.0;
    for (size_t i = 0; i < gc.size(); ++i) {
        dot_val += gc(i) * (this->rcsd(i) - r(i));
    }
    hc += dot_val;
    return {cut, shrunk, this->rcsd, this->_num_retries < 15};
}
