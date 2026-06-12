#include <ellalgo/arr.hpp>
#include <multiplierless/lowpass_oracle_q.hpp>

using Vec = std::valarray<double>;
using ParallelCut = std::pair<Arr, Vec>;

extern auto csd_quantize(double num, unsigned int) -> double;

extern auto inverse_spectral_fact(const Arr& r) -> Arr;
extern auto spectral_fact(const Arr& r) -> Arr;

/**
 * @brief Assess the optimization with quantized (CSD) coefficients.
 *
 * On the first call (retry=false), evaluates the lowpass oracle with the
 * continuous coefficients, performs spectral factorization, quantizes the
 * impulse response to CSD format, and computes the quantized autocorrelation.
 * On retry calls, evaluates the oracle with the previously computed quantized
 * coefficients and adjusts the cutting plane by the quantization error.
 *
 * @param[in]     r     The top-half autocorrelation coefficients (continuous).
 * @param[in,out] Spsq  On input, target stopband attenuation squared.
 *                      On output, achieved stopband value.
 * @param[in]     retry True if the previous cut had no effect and this is a retry.
 *
 * @return A tuple of (ParallelCut, is_optimal, rcsd, more_retries).
 *         - ParallelCut:  (gradient, objective) cutting-plane pair.
 *         - bool:         True if the optimal solution was found.
 *         - Arr:          The quantized autocorrelation coefficients used.
 *         - bool:         True if more retries should be attempted.
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
            hcsd(i) = csd_quantize(h(i), this->_nnz);
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
