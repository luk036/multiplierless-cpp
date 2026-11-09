#include <doctest/doctest.h>

#include <cmath>
#include <ellalgo/arr.hpp>
#include <multiplierless/spectral_fact.hpp>
#include <stdexcept>

TEST_CASE("spectral_fact round-trip (Python reference test)") {
    // Exact test case from Python's test_spectral_fact.py
    Arr h(7);
    h(0) = 0.76006445F;
    h(1) = 0.54101887F;
    h(2) = 0.42012073F;
    h(3) = 0.3157191F;
    h(4) = 0.10665804F;
    h(5) = 0.04326203F;
    h(6) = 0.01315678F;

    auto r = inverse_spectral_fact(h);
    auto h2 = spectral_fact_fft(r);

    CHECK(h.size() == h2.size());
    for (size_t i = 0; i < h.size(); ++i) {
        CHECK(h2(i) == doctest::Approx(h(i)).epsilon(0.05));
    }
}

TEST_CASE("spectral_fact rejects invalid autocorrelation") {
    CHECK_THROWS_AS(spectral_fact_fft(zeros(0)), std::invalid_argument);

    auto negative = zeros(4);
    negative(0) = -1.0;
    CHECK_THROWS_AS(spectral_fact_fft(negative), std::invalid_argument);

    auto non_finite = zeros(4);
    non_finite(0) = 1.0;
    non_finite(1) = std::nan("");
    CHECK_THROWS_AS(spectral_fact_fft(non_finite), std::invalid_argument);
}

TEST_CASE("spectral_fact_root rejects bad tolerance") {
    auto r = zeros(4);
    r(0) = 1.0;
    CHECK_THROWS_AS(spectral_fact_root(r, -1.0), std::invalid_argument);
}

TEST_CASE("spectral_fact rejects significantly negative spectrum") {
    auto r = zeros(3);
    r(0) = 0.01;
    r(1) = 0.5;
    r(2) = 0.3;
    CHECK_THROWS_AS(spectral_fact_fft(r), std::runtime_error);
}
