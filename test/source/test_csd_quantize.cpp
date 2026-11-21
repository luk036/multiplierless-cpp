#include <doctest/doctest.h>

#include <cmath>
#include <cstdlib>
#include <multiplierless/csd_quantize.hpp>
#include <multiplierless/verilog_utils.hpp>
#include <numbers>
#include <string>
#include <vector>

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
    auto q3_lo = csd_quantize(std::numbers::pi, 2);
    auto q3_hi = csd_quantize(std::numbers::pi, 8);
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

TEST_CASE("CSD string round-trips the quantized value") {
    for (auto const value : {0.0, 1.0, -1.0, 0.5, -0.5, 3.25, 0.9127555772777217}) {
        auto const parts = csd_quantize_parts(value, 7);
        auto const csd = format_csd(parts.digits);
        CHECK(csd.find('.') != std::string::npos);
        CHECK(csd_to_float(csd) == doctest::Approx(parts.value).epsilon(1e-15));
        CHECK(csd_quantize(value, 7) == doctest::Approx(parts.value));
    }
}

TEST_CASE("fixed-point coefficients share one scale") {
    std::vector<std::string> const csds{"0.000000+0+000+000+", "0.+"};
    auto const fixed = csd_coeffs_to_fixed(csds);
    CHECK(fixed.max_power == 17);
    for (size_t i = 0; i < fixed.digits.size(); ++i) {
        auto const& raw = fixed.digits[i];
        CHECK(raw.size() == static_cast<size_t>(fixed.max_power + 1));
        auto value = 0LL;
        for (size_t j = 0; j < raw.size(); ++j) {
            auto const power = static_cast<int>(raw.size() - 1 - j);
            if (raw[j] == '+') value += 1LL << power;
            if (raw[j] == '-') value -= 1LL << power;
        }
        auto const want = std::llround(csd_to_float(csds[i]) * std::pow(2.0, fixed.scale_power));
        CHECK(value == want);
    }
}
