/// @brief Quick smoke-test for spectral factorization.

#include <cstdio>
#include <ellalgo/arr.hpp>
#include <multiplierless/fftw_helper.hpp>

auto spectral_fact(const Arr& r) -> Arr;

/**
 * @brief Test spectral factorisation with a simple AR(1) autocorrelation.
 * @return 0 on success.
 */
int main() {
    // AR(1) with alpha=0.3, n=8
    Arr r(8);
    r(0) = 1.0;
    for (size_t k = 1; k < 8; ++k) r(k) = 0.3;
    printf("r =");
    for (size_t i = 0; i < 8; ++i) printf(" %.3f", r(i));
    printf("\n");

    auto h = spectral_fact(r);
    printf("h =");
    for (size_t i = 0; i < 8; ++i) printf(" %.6f", h(i));
    printf("\n");
    printf("OK\n");
    return 0;
}
