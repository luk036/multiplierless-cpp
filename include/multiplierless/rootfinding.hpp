#pragma once

// import numpy as np
#include <tuple>
#include <vector>

#include "matrix2.hpp"
#include "vector2.hpp"

using vec2 = numeric::vector2<double>;
using mat2 = numeric::matrix2<vec2>;

extern auto makeG(const vec2& vr, const vec2& vp) -> mat2;
extern auto makeadjoint(const vec2& vr, const vec2& vp) -> mat2;
extern void suppress(const vec2& vA, vec2& vA1, const vec2& vr, const vec2& vrj);
extern auto check_newton(const vec2& vA, const vec2& vA1, const vec2& vr) -> vec2;
extern auto horner_eval(std::vector<double>& pb, std::size_t n, const double& r) -> double;
extern auto horner(std::vector<double>& pb, std::size_t n, const vec2& vr) -> vec2;

class Options {
  public:
    unsigned int max_iter = 2000U;
    double tol = 1e-14;
};

extern auto initial_guess(const std::vector<double>& pa) -> std::vector<vec2>;
extern auto pbairstow_even(const std::vector<double>& pa, std::vector<vec2>& vrs,
                           const Options& options) -> std::tuple<unsigned int, bool>;
