#define DOCTEST_CONFIG_USE_STD_HEADERS
#include <doctest/doctest.h>

#include <ellalgo/arr.hpp>
#include <multiplierless/lowpass_oracle.hpp>

TEST_CASE("LowpassOracle assess_optim x(0) < 0 branch") {
    auto Fdc = filter_design_construct(16);
    auto Spsq = Fdc.Spsq;
    auto omega = LowpassOracle(std::move(Fdc));

    auto x = zeros(16);
    x(0) = -1.0;  // Trigger the x(0) < 0 branch

    auto [cut, feasible] = omega.assess_optim(x, Spsq);
    CHECK(!feasible);
    auto& [g, f] = cut;
    CHECK(g(0) == doctest::Approx(-1.0));
    CHECK(f.size() == 1);
    CHECK(f[0] == doctest::Approx(1.0));
}

TEST_CASE("LowpassOracle assess_optim initial feasible") {
    // With r(0) > 0, the oracle should proceed past the x(0) < 0 check
    auto Fdc = filter_design_construct(16);
    auto Spsq = Fdc.Spsq;
    auto omega = LowpassOracle(std::move(Fdc));

    auto x = zeros(16);
    x(0) = 1.0;  // Positive first coefficient

    // Should not hit the x(0) < 0 branch
    auto [cut, feasible] = omega.assess_optim(x, Spsq);
    // Should get some cut or feasible result without crashing
    CHECK(!feasible);  // unlikely to be optimal on first call
    CHECK(cut.first.size() > 0);
    CHECK(cut.second.size() > 0);
}

TEST_CASE("LowpassOracle assess_optim passband constraint check") {
    auto Fdc = filter_design_construct(8);
    auto Spsq = Fdc.Spsq;
    auto omega = LowpassOracle(std::move(Fdc));

    // A decent initial guess for autocorrelation coefficients
    auto x = zeros(8);
    for (size_t i = 0; i < 8; ++i) x(i) = 1.0 / static_cast<double>(i + 1);

    x(0) = 5.0;
    auto [cut, feasible] = omega.assess_optim(x, Spsq);
    CHECK(cut.first.size() > 0);
    CHECK(cut.second.size() > 0);
}

TEST_CASE("LowpassOracle operator() convenience") {
    auto Fdc = filter_design_construct(8);
    auto Spsq = Fdc.Spsq;
    auto omega = LowpassOracle(std::move(Fdc));

    auto x = zeros(8);
    x(0) = -2.0;

    auto [cut, feasible] = omega(x, Spsq);
    CHECK(!feasible);
    CHECK(cut.first(0) == doctest::Approx(-1.0));
    CHECK(cut.second[0] == doctest::Approx(2.0));
}
