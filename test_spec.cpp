#include <ellalgo/arr.hpp>
#include <multiplierless/fftw_helper.hpp>
#include <cstdio>

auto spectral_fact(const Arr& r) -> Arr;

int main() {
    // AR(1) with alpha=0.3, n=8
    Arr r(8);
    r(0) = 1.0;
    for (size_t k = 1; k < 8; ++k) r(k) = 0.3;
    printf("r ="); for (size_t i=0;i<8;++i) printf(" %.3f", r(i)); printf("\n");

    auto h = spectral_fact(r);
    printf("h ="); for (size_t i=0;i<8;++i) printf(" %.6f", h(i)); printf("\n");
    printf("OK\n");
    return 0;
}
