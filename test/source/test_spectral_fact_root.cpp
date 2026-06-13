#define DOCTEST_CONFIG_USE_STD_HEADERS
#include <doctest/doctest.h>

#include <cmath>
#include <ellalgo/arr.hpp>

extern auto spectral_fact_root(const Arr& r, double tolerance) -> Arr;

/**
 * @brief Fill r with a valid auto-correlation sequence: r(k) = alpha^k.
 * @param[in,out] r     Array to fill (size determines sequence length).
 * @param[in]     alpha  Decay factor in (0, 1).
 */
static void fill_acorr_decay(Arr& r, double alpha) {
    r(0) = 1.0;
    double ak = alpha;
    for (size_t k = 1; k < r.size(); ++k) {
        r(k) = ak;
        ak *= alpha;
    }
}

TEST_CASE("spectral_fact_root produces finite output") {
    auto n = static_cast<size_t>(16);
    auto r = Arr(n);
    fill_acorr_decay(r, 0.5);

    auto h = spectral_fact_root(r, 1e-10);
    CHECK(h.size() == n);
    for (size_t i = 0; i < n; ++i) CHECK(std::isfinite(h(i)));
}

TEST_CASE("spectral_fact_root energy matches r(0)") {
    auto n = static_cast<size_t>(16);
    auto r = Arr(n);
    fill_acorr_decay(r, 0.6);
    auto r0 = r(0);

    auto h = spectral_fact_root(r, 1e-10);
    double energy = 0.0;
    for (size_t i = 0; i < n; ++i) energy += h(i) * h(i);

    auto rel_err = std::fabs((energy - r0) / (r0 + 1e-10));
    CHECK(rel_err < 0.05);
}

TEST_CASE("spectral_fact_root with mild decay") {
    auto n = static_cast<size_t>(16);
    auto r = Arr(n);
    fill_acorr_decay(r, 0.8);

    auto h = spectral_fact_root(r, 1e-10);
    CHECK(h.size() == n);
    for (size_t i = 0; i < n; ++i) CHECK(std::isfinite(h(i)));

    double energy = 0.0;
    for (size_t i = 0; i < n; ++i) energy += h(i) * h(i);
    auto rel_err = std::fabs((energy - r(0)) / (r(0) + 1e-10));
    CHECK(rel_err < 0.05);
}

TEST_CASE("spectral_fact_root smaller order") {
    auto n = static_cast<size_t>(8);
    auto r = Arr(n);
    fill_acorr_decay(r, 0.5);

    auto h = spectral_fact_root(r, 1e-10);
    CHECK(h.size() == n);
    for (size_t i = 0; i < n; ++i) CHECK(std::isfinite(h(i)));
}
