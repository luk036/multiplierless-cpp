#define DOCTEST_CONFIG_USE_STD_HEADERS
#include <doctest/doctest.h>

#include <cmath>
#include <cstdlib>

extern auto csd_quantize(double num, unsigned int nnz) -> double;

TEST_CASE("csd_quantize zero input") {
    CHECK(csd_quantize(0.0, 0) == doctest::Approx(0.0));
    CHECK(csd_quantize(0.0, 1) == doctest::Approx(0.0));
    CHECK(csd_quantize(0.0, 10) == doctest::Approx(0.0));
}

TEST_CASE("csd_quantize positive values") {
    // A small positive number
    auto q1 = csd_quantize(0.5, 4);
    CHECK(q1 > 0.0);
    CHECK(std::fabs(q1 - 0.5) < 0.1);

    // A larger positive number
    auto q2 = csd_quantize(42.0, 4);
    CHECK(q2 > 0.0);
    CHECK(std::fabs(q2 - 42.0) < 5.0);

    // With more non-zero digits, should be more accurate
    auto q3_lo = csd_quantize(3.14159, 2);
    auto q3_hi = csd_quantize(3.14159, 8);
    CHECK(std::fabs(q3_hi - 3.14159) <= std::fabs(q3_lo - 3.14159) + 1e-12);
}

TEST_CASE("csd_quantize negative values") {
    auto q = csd_quantize(-0.5, 4);
    CHECK(q < 0.0);
    CHECK(std::fabs(q + 0.5) < 0.1);

    auto q2 = csd_quantize(-100.0, 6);
    CHECK(q2 < 0.0);
}

TEST_CASE("csd_quantize with nnz=0") {
    // With zero non-zero digits, result should be 0
    auto q = csd_quantize(100.0, 0);
    CHECK(q == doctest::Approx(0.0));
}

TEST_CASE("csd_quantize with nnz=1") {
    // With just 1 non-zero digit, should be a power-of-two approximation
    auto q = csd_quantize(1.5, 1);
    CHECK(q > 0.0);
    // Should be close to a power of 2 (1.0 or 2.0)
    CHECK((std::fabs(q - 1.0) < 0.01 || std::fabs(q - 2.0) < 0.01));
}

TEST_CASE("csd_quantize very small positive") {
    auto q = csd_quantize(1e-50, 4);
    CHECK(q == doctest::Approx(0.0));
}

TEST_CASE("csd_quantize large positive") {
    auto q = csd_quantize(1e6, 4);
    CHECK(q > 0.0);
    CHECK(q < 1.1e6);
}

TEST_CASE("csd_quantize symmetry") {
    // Positive and negative should be opposite
    auto pos = csd_quantize(1.234, 4);
    auto neg = csd_quantize(-1.234, 4);
    CHECK(pos + neg == doctest::Approx(0.0));
}
