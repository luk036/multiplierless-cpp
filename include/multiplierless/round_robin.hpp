// -*- coding: utf-8 -*-
#pragma once

/** @file round_robin.hpp
 *  @brief Round-robin index scanner for constraint cycling.
 */

#include <cstddef>

/**
 * @brief Round-robin counter over the half-open range [lo, hi).
 *
 * Each call to next() returns the next index in the cycle, wrapping back to
 * lo after reaching hi. Used to scan constraint rows in rotating order so
 * that no single constraint starves the cutting-plane method.
 */
class RoundRobin {
    std::size_t _lo;
    std::size_t _hi;
    std::size_t _cur;

  public:
    /// @brief Default constructor (empty range).
    RoundRobin() : _lo{0}, _hi{0}, _cur{0} {}

    /**
     * @brief Construct a round-robin counter over [lo, hi).
     *
     * The first call to next() returns lo.
     *
     * @param[in] lo Lower bound (inclusive).
     * @param[in] hi Upper bound (exclusive).
     */
    RoundRobin(std::size_t lo, std::size_t hi) : _lo{lo}, _hi{hi}, _cur{hi - 1} {}

    /**
     * @brief Return the next index in the cycle.
     *
     * @return The next index in [lo, hi), wrapping to lo after hi - 1.
     */
    auto next() -> std::size_t {
        if (++_cur == _hi) {
            _cur = _lo;
        }
        return _cur;
    }
};
