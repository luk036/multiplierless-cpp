#include <doctest/doctest.h>

#include <ellalgo/arr.hpp>

extern auto spectral_fact(const Arr& r) -> Arr;
extern auto inverse_spectral_fact(const Arr& h) -> Arr;

TEST_CASE("spectral_fact round-trip (Python reference test)") {
    // Exact test case from Python's test_spectral_fact.py
    Arr h(7);
    h(0) = 0.76006445f;
    h(1) = 0.54101887f;
    h(2) = 0.42012073f;
    h(3) = 0.3157191f;
    h(4) = 0.10665804f;
    h(5) = 0.04326203f;
    h(6) = 0.01315678f;

    auto r = inverse_spectral_fact(h);
    auto h2 = spectral_fact(r);

    CHECK(h.size() == h2.size());
    for (size_t i = 0; i < h.size(); ++i) {
        CHECK(h2(i) == doctest::Approx(h(i)).epsilon(1e-5));
    }
}
