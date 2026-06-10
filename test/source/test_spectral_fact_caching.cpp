#include <doctest/doctest.h>

#include <ellalgo/arr.hpp>

#include <cmath>

extern auto spectral_fact(const Arr& r) -> Arr;
extern auto spectral_fact_fft(const Arr& r) -> Arr;
extern auto inverse_spectral_fact(const Arr& h) -> Arr;

/// Fill r with a valid auto-correlation sequence: r(k) = alpha^k, alpha in (0,1)
static void fill_ar1_acorr(Arr& r, double alpha) {
    r(0) = 1.0;
    double ak = alpha;
    for (size_t k = 1; k < r.size(); ++k) {
        r(k) = ak;
        ak *= alpha;
    }
}

TEST_CASE("spectral_fact_fft caching with same N") {
    // First call builds the cache
    auto r1 = Arr(32);
    fill_ar1_acorr(r1, 0.5);

    auto h1 = spectral_fact_fft(r1);
    CHECK(h1.size() == 32);

    // Second call with same N hits the cache
    auto r2 = Arr(32);
    fill_ar1_acorr(r2, 0.6);

    auto h2 = spectral_fact_fft(r2);
    CHECK(h2.size() == 32);

    // Both should produce finite results
    for (size_t i = 0; i < 32; ++i) {
        CHECK(std::isfinite(h1(i)));
        CHECK(std::isfinite(h2(i)));
    }
}

TEST_CASE("spectral_fact_fft different orders swap cache") {
    auto r16 = Arr(16);
    fill_ar1_acorr(r16, 0.5);
    auto h16 = spectral_fact_fft(r16);
    CHECK(h16.size() == 16);

    auto r32 = Arr(32);
    fill_ar1_acorr(r32, 0.5);
    auto h32 = spectral_fact_fft(r32);
    CHECK(h32.size() == 32);

    // Calling 16 again re-builds cache
    auto h16b = spectral_fact_fft(r16);
    CHECK(h16b.size() == 16);
}

TEST_CASE("spectral_fact_fft round-trip via inverse") {
    auto n = static_cast<size_t>(16);
    auto r = Arr(n);
    fill_ar1_acorr(r, 0.5);
    auto r0 = r(0);

    auto h = spectral_fact_fft(r);
    auto r_back = inverse_spectral_fact(h);

    CHECK(r_back.size() == n);
    double energy = 0.0;
    for (size_t i = 0; i < n; ++i) energy += h(i) * h(i);
    CHECK(std::fabs(energy - r0) / (r0 + 1e-10) < 0.01);
}

TEST_CASE("spectral_fact convenience wrapper") {
    auto r = Arr(16);
    fill_ar1_acorr(r, 0.5);

    auto h = spectral_fact(r);
    CHECK(h.size() == 16);
    for (size_t i = 0; i < 16; ++i) CHECK(std::isfinite(h(i)));
}

TEST_CASE("inverse_spectral_fact with zero input") {
    auto h = zeros(8);
    auto r = inverse_spectral_fact(h);
    CHECK(r.size() == 8);
    for (size_t i = 0; i < 8; ++i) CHECK(r(i) == doctest::Approx(0.0));
}

TEST_CASE("inverse_spectral_fact with single impulse") {
    auto h = Arr(8);
    h(0) = 0.0;
    h(1) = 1.0;
    auto r = inverse_spectral_fact(h);
    CHECK(r.size() == 8);
    CHECK(r(0) == doctest::Approx(1.0));
}
