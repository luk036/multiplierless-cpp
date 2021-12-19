// import numpy as np
#include <cmath>  // import pow, cos, sqrt
#include <multiplierless/rootfinding.hpp>

// using vec2 = numeric::vector2<double>;
// using mat2 = numeric::matrix2<vec2>;

auto makeG(const vec2& vr, const vec2& vp) -> mat2 {
    const auto& [r, q] = vr;
    const auto& [p, s] = vp;
    return mat2{vec2{p * r + s, p}, vec2{p * q, s}};
}

auto makeadjoint(const vec2& vr, const vec2& vp) -> mat2 {
    const auto& [r, q] = vr;
    const auto& [p, s] = vp;
    return mat2{vec2{s, -p}, vec2{-p * q, p * r + s}};
}

void suppress(const vec2& vA, vec2& vA1, const vec2& vr, const vec2& vrj) {
    auto vp = vr - vrj;
    auto mp = makeadjoint(vrj, vp);  // 2 mul's
    vA1 -= mp.mdot(vA) / mp.det();  // 6 mul's + 2 div's
    return;
}

auto check_newton(const vec2& vA, const vec2& vA1, const vec2& vr) -> vec2 {
    auto mA1 = makeadjoint(vr, vA1); // 2 mul's
    return mA1.mdot(vA) / mA1.det(); // 6 mul's + 2 div's
}


auto horner_eval(std::vector<double>& pb, size_t n, const double& r) -> double {
    for (auto i = 0U; i != n; ++i) {
        pb[i + 1] += pb[i] * r;
    }
    return pb[n];
}

auto horner(std::vector<double>& pb, size_t n, const vec2& vr) -> vec2 {
    const auto& [r, q] = vr;
    pb[1] += pb[0] * r;
    for (auto i = 2U; i != n; ++i) {
        pb[i] += pb[i - 2] * q + pb[i - 1] * r;
    }
    pb[n] += pb[n - 2] * q;
    return vec2{pb[n - 1], pb[n]};
}

auto initial_guess(const std::vector<double>& pa) -> std::vector<vec2> {
    static const auto PI = std::acos(-1.);

    auto N = pa.size() - 1;
    auto M = N / 2;
    auto c = -pa[1]/(N*pa[0]);
    auto pb = pa;
    auto Pc = horner_eval(pb, N, c); // ???
    auto re = std::pow(std::abs(Pc), 1./N);
    auto k = 2 * PI / N;
    auto m = c * c + re * re;
    auto vr0s = std::vector<vec2>{};
    for (auto i = 1U; i != M + 1; ++i) {
        auto r0 = 2 * (c + re * std::cos(k * i));
        auto q0 = m + r0;
        vr0s.emplace_back(vec2{r0, q0});
    }
    return vr0s;
}

auto pbairstow_even(const std::vector<double>& pa, std::vector<vec2>& vrs, const Options& options = Options()) -> std::tuple<unsigned int, bool> {
    auto N = pa.size() - 1; // degree, assume even
    auto M = N / 2;
    auto found = false;
    auto niter = 0U;
    for (; niter != options.max_iter; ++niter) {
        auto tol = 0.;
        for (auto i = 0U; i != M; ++i) {
            auto pb = pa;
            // auto n = pa.size() - 1;
            auto vA = horner(pb, N, vrs[i]);
            const auto& [A, B] = vA;
            auto toli = std::abs(A) + std::abs(B);
            if (toli < options.tol) {
                continue;
            }
            tol = std::max(tol, toli);
            auto vA1 = horner(pb, N - 2, vrs[i]);
            for (auto j = 0U; j != M; ++j) { // exclude i
                if (j == i) {
                    continue;
                }
                auto vp = vrs[i] - vrs[j];
                auto mp = makeadjoint(vrs[j], vp);  // 2 mul's
                vA1 -= mp.mdot(vA) / mp.det();  // 6 mul's + 2 div's
                // vA1 = suppress(vA, vA1, vrs[i], vrs[j]);
            }
            auto mA1 = makeadjoint(vrs[i], vA1); // 2 mul's
            vrs[i] -= mA1.mdot(vA) / mA1.det(); // Gauss-Seidel fashion
        }
        // fmt::print("tol: {}\n", tol);
        if (tol < options.tol) {
            found = true;
            break;
        }
    }
    return {++niter, found};
}

// auto find_rootq(const vec2& r) {
//     auto hb = b / 2.;
//     auto d = hb * hb - c;
//     if (d < 0.) {
//         auto x1 = -hb + (sqrt(-d) if (hb < 0. else -sqrt(-d))*1j;
//     }
//     else {
//         auto x1 = -hb + (sqrt(d) if (hb < 0. else -sqrt(d));
//     }
//     auto x2 = c / x1;
//     return x1, x2;
// }

