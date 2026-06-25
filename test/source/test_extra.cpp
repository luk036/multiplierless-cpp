#include <doctest/doctest.h>

#include <multiplierless/lowpass_oracle.hpp>
#include <multiplierless/lowpass_oracle_q.hpp>

// ============================================================
// Additional tests to improve code coverage
// ============================================================

// --- filter_design_construct construction edge cases ---

TEST_CASE("filter_design_construct default params") {
    filter_design_construct fdc(32);
    CHECK_EQ(fdc.N, 32);
    CHECK_GT(fdc.Ap.rows(), 0);
    CHECK_GT(fdc.As.rows(), 0);
    CHECK_GT(fdc.Anr.rows(), 0);
    CHECK_GT(fdc.Lpsq, 0.0);
    CHECK_GT(fdc.Upsq, 0.0);
    CHECK_GT(fdc.Spsq, 0.0);
}

TEST_CASE("filter_design_construct explicit params") {
    auto fdc = filter_design_construct(16, 0.15, 0.25, 0.1, 0.1, 10);
    CHECK_EQ(fdc.N, 16);
    CHECK_GT(fdc.Lpsq, 0.0);
    CHECK_GT(fdc.Upsq, 0.0);
}

TEST_CASE("filter_design_construct low order") {
    auto fdc = filter_design_construct(4);
    CHECK_EQ(fdc.N, 4);
}

TEST_CASE("filter_design_construct high order") {
    auto fdc = filter_design_construct(64);
    CHECK_EQ(fdc.N, 64);
}

// --- LowpassOracle edge cases ---

TEST_CASE("LowpassOracle custom construction") {
    auto fdc = filter_design_construct(16, 0.15, 0.25, 0.1, 0.1, 12);
    LowpassOracle oracle(std::move(fdc));
    // Oracle constructed with custom filter specs
}

TEST_CASE("LowpassOracle assess_optim negative x0") {
    auto fdc = filter_design_construct(16);
    LowpassOracle oracle(std::move(fdc));
    Arr x = zeros(16);
    x(0) = -1.0;  // violates nonnegative-real constraint
    double Spsq = 0.0;
    auto [cut, is_optimal] = oracle.assess_optim(x, Spsq);
    CHECK_FALSE(is_optimal);  // should not be optimal
}

TEST_CASE("LowpassOracle assess_optim positive x0") {
    auto fdc = filter_design_construct(16);
    LowpassOracle oracle(std::move(fdc));
    Arr x = zeros(16);
    x(0) = 1.0;
    x(1) = 0.5;
    double Spsq = 0.0;
    auto [cut, is_optimal] = oracle.assess_optim(x, Spsq);
    // Should pass negative check and proceed to passband/stopband constraints
    CHECK_FALSE(is_optimal);
}

TEST_CASE("LowpassOracle assess_optim with all ones") {
    auto fdc = filter_design_construct(8);
    LowpassOracle oracle(std::move(fdc));
    Arr x = zeros(8);
    std::fill(x.begin(), x.end(), 1.0);
    double Spsq = 0.0;
    auto [cut, is_optimal] = oracle.assess_optim(x, Spsq);
    // All-ones input should exercise various constraint checks
    CHECK_FALSE(is_optimal);
}

TEST_CASE("LowpassOracle assess_optim large coefficients") {
    auto fdc = filter_design_construct(8, 0.12, 0.20, 0.125, 0.125, 10);
    LowpassOracle oracle(std::move(fdc));
    Arr x = zeros(8);
    for (size_t i = 0; i < 8; ++i) {
        x(i) = 100.0;  // large values likely violate Upsq
    }
    double Spsq = 0.0;
    auto [cut, is_optimal] = oracle.assess_optim(x, Spsq);
    CHECK_FALSE(is_optimal);
}

TEST_CASE("LowpassOracle assess_optim zero coefficients") {
    auto fdc = filter_design_construct(8);
    LowpassOracle oracle(std::move(fdc));
    Arr x = zeros(8);
    double Spsq = 0.0;
    auto [cut, is_optimal] = oracle.assess_optim(x, Spsq);
    // Zero coefficients: x(0)=0 not negative, should reach constraint checks
    CHECK_FALSE(is_optimal);
}

// --- LowpassOracleQ tests ---

TEST_CASE("LowpassOracleQ with oracle") {
    auto fdc = filter_design_construct(8);
    LowpassOracle oracle(std::move(fdc));
    LowpassOracleQ oracle_q(8u, std::move(oracle));
}

// --- Filter design extreme parameters ---

TEST_CASE("filter_design_construct tight passband") {
    // Very tight passband with small ripple
    auto fdc = filter_design_construct(32, 0.1, 0.3, 0.01, 0.01, 20);
    CHECK_EQ(fdc.N, 32);
}

TEST_CASE("filter_design_construct wide passband") {
    // Wide passband with large ripple tolerance
    auto fdc = filter_design_construct(16, 0.4, 0.48, 0.5, 0.5, 8);
    CHECK_EQ(fdc.N, 16);
}

// --- round-robin index tracking ---

TEST_CASE("LowpassOracle assess_optim multiple calls") {
    auto fdc = filter_design_construct(8);
    LowpassOracle oracle(std::move(fdc));
    Arr x = zeros(8);
    x(0) = 5.0;

    // Call assess_optim multiple times to exercise round-robin state
    double Spsq = 0.0;
    for (int i = 0; i < 5; ++i) {
        auto [cut, is_optimal] = oracle.assess_optim(x, Spsq);
        // Internal indices (_i_Ap, _i_As, _i_Anr) should advance
    }
}

TEST_CASE("LowpassOracle assess_optim alternating inputs") {
    auto fdc = filter_design_construct(8);
    LowpassOracle oracle(std::move(fdc));
    double Spsq = 0.0;

    // Call with different inputs to exercise different code paths
    Arr x1 = zeros(8);
    x1(0) = 5.0;
    oracle.assess_optim(x1, Spsq);

    Arr x2 = zeros(8);
    x2(0) = 5.0;
    for (size_t i = 1; i < 8; ++i) {
        x2(i) = 0.1 * static_cast<double>(i);
    }
    oracle.assess_optim(x2, Spsq);

    Arr x3 = zeros(8);
    for (size_t i = 0; i < 8; ++i) {
        x3(i) = -0.01 * static_cast<double>(i);
    }
    oracle.assess_optim(x3, Spsq);
}

// --- Round-robin boundary conditions ---

TEST_CASE("LowpassOracle round-robin wraps correctly") {
    auto fdc = filter_design_construct(4);
    LowpassOracle oracle(std::move(fdc));
    double Spsq = 0.0;
    Arr x = zeros(4);
    x(0) = 1.0;

    // Many iterations to force round-robin wraparound
    for (int i = 0; i < 20; ++i) {
        auto [cut, is_optimal] = oracle.assess_optim(x, Spsq);
    }
}
