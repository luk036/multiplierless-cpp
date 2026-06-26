#pragma once

/** @file lowpass_oracle_q.hpp
 *  @brief Oracle for quantized FIR lowpass filter design with CSD representation.
 */

#include <ellalgo/arr.hpp>
#include <multiplierless/lowpass_oracle.hpp>

/*!
 * @brief Oracle for quantized FIR lowpass filter design.
 *
 * This class wraps LowpassOracle and adds support for quantized (multiplierless)
 * filter coefficients using Canonical Signed Digit (CSD) representation.
 * It converts the filter coefficients to CSD format and performs optimization
 * with the quantized values while tracking the number of retries.
 */
class LowpassOracleQ {
    using Vec = std::valarray<double>;
    using ParallelCut = std::pair<Arr, Vec>;

    Arr rcsd{};                     ///< Quantized filter coefficients in CSD representation
    unsigned int _nnz;              ///< Maximum number of non-zero CSD digits
    LowpassOracle _lowpass;         ///< The underlying lowpass oracle for optimization
    unsigned int _num_retries = 0;  ///< Number of retry attempts made

  public:
    /*!
     * @brief Construct a new LowpassOracleQ object.
     *
     * @param[in] nnz Maximum number of non-zero digits allowed in CSD representation.
     * @param[in] lowpass The LowpassOracle instance to wrap for optimization.
     */
    LowpassOracleQ(unsigned int nnz, LowpassOracle&& lowpass)
        : _nnz(nnz), _lowpass(std::move(lowpass)) {}

    /*!
     * @brief Assess optimization with quantized filter coefficients.
     *
     * This function evaluates the optimization problem using quantized coefficients
     * represented in Canonical Signed Digit (CSD) format. It converts the input
     * coefficients to CSD, performs the optimization assessment, and returns the
     * cutting plane information along with the quantized coefficients.
     *
     * @param[in] r The filter coefficients (autocorrelation coefficients).
     *
     * @f[
     *     \begin{array}{ll}
     *     \text{minimize} & \max_{\omega \in \Omega_s} R(\omega) \\
     *     \text{subject to} & L^2 \le R(\omega) \le U^2,\; \omega \in \Omega_p \\
     *                       & R(\omega) \ge 0,\; r \in \mathcal{C}
     *     \end{array}
     * @f]
     * where \f$R(\omega) = r_0 + 2\sum_{k=1}^{N-1} r_k \cos(k\omega)\f$ and
     * \f$\mathcal{C}\f$ denotes the CSD quantization constraint.
     *
     * @param[in,out] Spsq On input, the target stopband attenuation squared.
     *                      On output, the achieved maximum stopband value.
     * @param[in] retry Whether this is a retry after the previous cut had no effect.
     *
     * @return A tuple containing:
     *         - ParallelCut: Pair of gradient and objective function values
     *         - bool: True if optimal solution found
     *         - Arr: The quantized filter coefficients used
     *         - bool: True if more retries should be attempted
     */
    auto assess_optim_q(const Arr& r, double& Spsq, bool retry)
        -> std::tuple<ParallelCut, bool, Arr, bool>;

    /*!
     * @brief Operator function for optimization assessment.
     *
     * Convenience function that forwards to assess_optim_q.
     *
     * @param[in] r The filter coefficients (autocorrelation coefficients).
     * @param[in,out] Spsq On input, the target stopband attenuation squared.
     *                      On output, the achieved maximum stopband value.
     * @param[in] retry Whether this is a retry after the previous cut had no effect.
     *
     * @return A tuple containing:
     *         - ParallelCut: Pair of gradient and objective function values
     *         - bool: True if optimal solution found
     *         - Arr: The quantized filter coefficients used
     *         - bool: True if more retries should be attempted
     */
    auto operator()(const Arr& r, double& Spsq, bool retry)
        -> std::tuple<ParallelCut, bool, Arr, bool> {
        return this->assess_optim_q(r, Spsq, retry);
    }
};
