#define DOCTEST_CONFIG_USE_STD_HEADERS
#include <doctest/doctest.h>

#include <ellalgo/arr.hpp>
#include <multiplierless/lowpass_oracle.hpp>

TEST_CASE("filter_design_construct default N=32") {
    auto fdc = filter_design_construct(32);
    CHECK(fdc.N == 32);
    CHECK(fdc.Lpsq > 0.0);
    CHECK(fdc.Upsq > 0.0);
    CHECK(fdc.Spsq > 0.0);
    CHECK(fdc.Ap.rows() > 0);
    CHECK(fdc.Ap.cols() == static_cast<size_t>(fdc.N));
    CHECK(fdc.As.rows() > 0);
    CHECK(fdc.Anr.rows() > 0);
}

TEST_CASE("filter_design_construct custom parameters") {
    // Custom passband 0-0.15π, stopband 0.25-π
    auto fdc = filter_design_construct(32, 0.15, 0.25, 0.1, 0.1, 10);
    CHECK(fdc.N == 32);
    CHECK(fdc.Lpsq > 0.0);
    CHECK(fdc.Upsq > fdc.Lpsq);
    CHECK(fdc.Spsq > 0.0);
    CHECK(fdc.Ap.rows() > 0);
    CHECK(fdc.As.rows() > 0);
    CHECK(fdc.Anr.rows() > 0);
}

TEST_CASE("filter_design_construct smaller order") {
    auto fdc = filter_design_construct(8);
    CHECK(fdc.N == 8);
    CHECK(fdc.Ap.cols() == 8);
}

TEST_CASE("filter_design_construct different passband stopband") {
    // Narrow passband, wide stopband
    auto fdc = filter_design_construct(16, 0.08, 0.30, 0.05, 0.05, 12);
    CHECK(fdc.N == 16);
    CHECK(fdc.Lpsq > 0.0);
    CHECK(fdc.Upsq > 0.0);
    CHECK(fdc.Spsq > 0.0);

    // Wide passband, narrow stopband
    auto fdc2 = filter_design_construct(16, 0.30, 0.40, 0.2, 0.01, 12);
    CHECK(fdc2.N == 16);
    CHECK(fdc2.Lpsq > 0.0);
    CHECK(fdc2.Upsq > 0.0);
    CHECK(fdc2.Spsq > 0.0);
}

TEST_CASE("filter_design_construct tight ripple specs") {
    // Very tight ripple specs
    auto fdc = filter_design_construct(64, 0.2, 0.24, 0.01, 0.01, 20);
    CHECK(fdc.N == 64);
    CHECK(fdc.Lpsq > 0.0);
    CHECK(fdc.Upsq > fdc.Lpsq);
    CHECK(fdc.Ap.rows() > 0);
    CHECK(fdc.As.rows() > 0);
}
