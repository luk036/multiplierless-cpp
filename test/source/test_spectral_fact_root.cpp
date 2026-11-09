#include <doctest/doctest.h>

#include <cmath>
#include <ellalgo/arr.hpp>
#include <multiplierless/spectral_fact.hpp>

#ifdef NDEBUG

TEST_CASE("spectral_fact_root preserves energy") {
    Arr h(7);
    h(0) = 0.76006445;
    h(1) = 0.54101887;
    h(2) = 0.42012073;
    h(3) = 0.3157191;
    h(4) = 0.10665804;
    h(5) = 0.04326203;
    h(6) = 0.01315678;

    auto const r = inverse_spectral_fact(h);
    auto const hr = spectral_fact_root(r, 1e-8);

    CHECK(hr.size() == h.size());
    auto energy = 0.0;
    for (size_t i = 0; i < hr.size(); ++i) energy += hr(i) * hr(i);
    CHECK(energy == doctest::Approx(r(0)).epsilon(0.05));
}

#endif
