// import numpy as np
// -*- coding: utf-8 -*-
#include <doctest/doctest.h>  // for ResultBuilder, CHECK, TEST...
#include <fmt/format.h>       // for print

#include <multiplierless/rootfinding.hpp>  // for horner, initial_guess, pba...
#include <vector>                          // for vector

#include "multiplierless/vector2.hpp"  // for vector2

TEST_CASE("test root-finding") {
    // auto vA = vec2{0.1, 1.2};
    // auto vA1 = vec2{2.3, 3.4};
    // auto vr = vec2{4.5, 5.6};
    // auto vrj = vec2{6.7, 7.8};
    // auto vA1 = suppress(vA, vA1, vr, vrj);
    // fmt::print(check_newton(vA, vA1, vr));
    auto h = std::vector<double>{5.0, 2.0, 9.0, 6.0, 2.0};
    auto vrs = initial_guess(h);
    // fmt::print(vrs);
    fmt::print("vrs[1]: {}, {}\n", vrs[1].x(), vrs[1].y());
    auto pb = h;
    auto N = pb.size() - 1;
    auto vAh = horner(pb, N, vrs[1]);
    fmt::print("{}, {}\n", vAh.x(), vAh.y());
    // fmt::print(pb);
    auto vA1h = horner(pb, N - 2, vrs[1]);
    fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

    auto [niter, found] = pbairstow_even(h, vrs, Options());
    fmt::print("{}, {}\n", niter, found);

    CHECK(niter <= 8);
    // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
}

TEST_CASE("test root-finding 2") {
    // auto vA = vec2{0.1, 1.2};
    // auto vA1 = vec2{2.3, 3.4};
    // auto vr = vec2{4.5, 5.6};
    // auto vrj = vec2{6.7, 7.8};
    // auto vA1 = suppress(vA, vA1, vr, vrj);
    // fmt::print(check_newton(vA, vA1, vr));
    auto h = std::vector<double>{10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0};
    auto vrs = initial_guess(h);
    // fmt::print(vrs);
    fmt::print("vrs[1]: {}, {}\n", vrs[1].x(), vrs[1].y());
    auto pb = h;
    auto N = pb.size() - 1;
    auto vAh = horner(pb, N, vrs[1]);
    fmt::print("{}, {}\n", vAh.x(), vAh.y());
    // fmt::print(pb);
    auto vA1h = horner(pb, N - 2, vrs[1]);
    fmt::print("{}, {}\n", vA1h.x(), vA1h.y());

    auto [niter, found] = pbairstow_even(h, vrs, Options());
    fmt::print("{}, {}\n", niter, found);

    CHECK(niter <= 13);
    // fmt::print([find_rootq(-r[0], -r[1]) for r : vrs]);
}
