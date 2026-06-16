#pragma once

/** @file logger.hpp
 *  @brief Logging utilities for the multiplierless library using spdlog.
 */

#include <string>

namespace multiplierless {

    /**
     * @brief Log a message using spdlog
     *
     * This function provides a simple wrapper around spdlog for logging messages.
     * It creates a file logger that writes to "multiplierless.log" and logs at the info level.
     *
     * @param message The message to log
     */
    void log_with_spdlog(const std::string& message);

}  // namespace multiplierless