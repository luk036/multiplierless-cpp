#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include <ellalgo/arr.hpp>
#include <ginger/aberth.hpp>
#include <ginger/config.hpp>
#include <multiplierless/fftw_helper.hpp>

#ifndef M_PI
constexpr double M_PI = 3.14159265358979323846264338327950288;
#endif

auto spectral_fact_fft(const Arr& r) -> Arr;

auto spectral_fact_root(const Arr& r, double tolerance) -> Arr {
    const auto n = r.size();
    const auto deg = 2 * n - 2;
    std::vector<double> coeffs(deg + 1, 0.0);
    coeffs[0] = r(n - 1);
    for (size_t i = 0; i < n - 1; ++i)
        coeffs[i + 1] = 2.0 * r(n - 2 - i);
    for (size_t i = 0; i < n - 2; ++i)
        coeffs[deg - i - 1] = 2.0 * r(n - 2 - i);
    coeffs[n - 1] = 2.0 * r(0);
    coeffs[deg] = r(n - 1);
    std::reverse(coeffs.begin(), coeffs.end());

    auto zs = initial_aberth_autocorr(coeffs);
    Options opts;
    opts.tolerance = tolerance;
    opts.max_iters = 500;
    aberth_autocorr(coeffs, zs, opts);

    std::vector<std::complex<double>> inside;
    for (auto& z : zs) {
        if (std::abs(z) < 1.0) inside.push_back(z);
        else inside.push_back(1.0 / z);
    }

    auto hc = poly_from_roots(inside);
    double eh = 0.0;
    for (auto c : hc) eh += c * c;
    const auto norm = std::sqrt(r(0) / eh);
    for (auto& c : hc) c *= norm;

    Arr h(n);
    for (size_t i = 0; i < n && i < hc.size(); ++i) h(i) = hc[i];
    return h;
}

auto spectral_fact(const Arr& r) -> Arr { return spectral_fact_fft(r); }

auto spectral_fact_fft(const Arr& r) -> Arr {
    const auto n = int(r.size());
    const auto mult_factor = 100;
    const auto m = mult_factor * n;

    static int cached_n = 0;
    static Arr cached_A;
    if (n != cached_n) {
        const auto step_w = 2.0 * M_PI / static_cast<double>(m);
        Arr w(m);
        for (size_t i = 0; i < size_t(m); ++i) w(i) = static_cast<double>(i) * step_w;
        auto cols = arange(1.0, double(n));
        Arr An = 2.0 * cos(outer(w, cols));
        cached_A = concatenate(ones(m, 1), An, 1);
        cached_n = n;
    }
    const auto& A = cached_A;
    Arr R = dot(A, r);

    auto min_val = *std::min_element(R.begin(), R.end());
    if (min_val <= 0) {
        for (size_t i = 0; i < R.size(); ++i)
            if (R(i) <= 0) R(i) = 1e-10;
    }

    Arr alpha = 0.5 * log(abs(R));
    auto alphatmp = fft(cast_to_complex(alpha));
    auto ind = size_t(m) / 2;
    for (auto i = ind; i < m; ++i) alphatmp[i] = -alphatmp[i];
    alphatmp[0] = {0.0, 0.0};
    alphatmp[ind] = {0.0, 0.0};

    const std::complex<double> j_{0, 1};
    Arr phi = ifft(j_ * alphatmp);
    auto alpha1 = view(alpha, Range(0, m, size_t(mult_factor)));
    auto phi1 = view(phi, Range(0, m, size_t(mult_factor)));
    return ifft(exp(cast_to_complex(alpha1) + j_ * cast_to_complex(phi1)));
}

auto inverse_spectral_fact(const Arr& h) -> Arr {
    auto n = h.size();
    Arr r(n);
    for (size_t t = 0; t < n; ++t) {
        double sum = 0.0;
        for (size_t i = 0; i < n - t; ++i) sum += h(i + t) * h(i);
        r(t) = sum;
    }
    return r;
}
