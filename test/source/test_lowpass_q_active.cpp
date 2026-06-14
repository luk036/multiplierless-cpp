#include <doctest/doctest.h>

#include <ellalgo/arr.hpp>
#include <multiplierless/lowpass_oracle_q.hpp>

class LowpassOracle;

extern auto create_lowpass_case(int N) -> std::tuple<LowpassOracle, double>;
extern auto csd_quantize(double num, unsigned int nnz) -> double;

TEST_CASE("LowpassOracleQ construction") {
    auto [omega, Spsq] = create_lowpass_case(16);
    auto omega_q = LowpassOracleQ(4, std::move(omega));
    CHECK(true);  // construction succeeded without exception
}

TEST_CASE("LowpassOracleQ assess_optim_q first call") {
    auto [omega, Spsq] = create_lowpass_case(16);
    auto omega_q = LowpassOracleQ(7, std::move(omega));

    auto r = zeros(16);
    r(0) = 1.0;
    // Seed with decaying values to make it a plausible autocorrelation
    for (size_t i = 1; i < 16; ++i) r(i) = 0.5 / static_cast<double>(i + 1);

    auto Spsq_val = Spsq;
    auto retry = false;

    auto [cut, shrunk, r_q, more_alt] = omega_q.assess_optim_q(r, Spsq_val, retry);
    CHECK(cut.first.size() > 0);
    CHECK(cut.second.size() > 0);
}

TEST_CASE("LowpassOracleQ assess_optim_q then retry") {
    auto [omega, Spsq] = create_lowpass_case(16);
    auto omega_q = LowpassOracleQ(7, std::move(omega));

    auto r = zeros(16);
    r(0) = 1.0;
    for (size_t i = 1; i < 16; ++i) r(i) = 0.5 / static_cast<double>(i + 1);

    // First call (non-retry) populates rcsd
    auto Spsq_val = Spsq;
    auto [cut1, shrunk1, r_q1, more_alt1] = omega_q.assess_optim_q(r, Spsq_val, false);
    CHECK(cut1.first.size() > 0);
    CHECK(cut1.second.size() > 0);

    // Second call with retry=true now has rcsd available
    if (shrunk1) {
        auto [cut2, shrunk2, r_q2, more_alt2] = omega_q.assess_optim_q(r, Spsq_val, true);
        CHECK(cut2.first.size() > 0);
        CHECK(cut2.second.size() > 0);
    }
}

TEST_CASE("LowpassOracleQ operator() convenience") {
    auto [omega, Spsq] = create_lowpass_case(16);
    auto omega_q = LowpassOracleQ(7, std::move(omega));

    auto r = zeros(16);
    r(0) = 1.0;
    for (size_t i = 1; i < 16; ++i) r(i) = 0.5 / static_cast<double>(i + 1);

    auto Spsq_val = Spsq;
    auto [cut, shrunk, r_q, more_alt] = omega_q(r, Spsq_val, false);
    CHECK(cut.first.size() > 0);
    CHECK(cut.second.size() > 0);
}
