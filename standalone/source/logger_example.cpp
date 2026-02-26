#include <multiplierless/logger.hpp>
#include <iostream>

/**
 * @brief Example program demonstrating spdlogger integration
 *
 * This program shows how to use the multiplierless::log_with_spdlog() function
 * to log messages to a file.
 */
int main() {
    std::cout << "Multiplierless Spdlogger Example" << std::endl;
    std::cout << "=================================" << std::endl;

    // Example 1: Basic logging
    std::cout << "\n[Example 1] Basic logging" << std::endl;
    multiplierless::log_with_spdlog("Application started");
    std::cout << "  -> Logged: Application started" << std::endl;

    // Example 2: Logging operational messages
    std::cout << "\n[Example 2] Logging operational messages" << std::endl;
    multiplierless::log_with_spdlog("Processing data...");
    std::cout << "  -> Logged: Processing data..." << std::endl;

    multiplierless::log_with_spdlog("Data processing completed successfully");
    std::cout << "  -> Logged: Data processing completed successfully" << std::endl;

    // Example 3: Logging status updates
    std::cout << "\n[Example 3] Logging status updates" << std::endl;
    multiplierless::log_with_spdlog("System initialized");
    std::cout << "  -> Logged: System initialized" << std::endl;

    multiplierless::log_with_spdlog("Ready to accept connections");
    std::cout << "  -> Logged: Ready to accept connections" << std::endl;

    std::cout << "\n=================================" << std::endl;
    std::cout << "Check 'multiplierless.log' for logged messages" << std::endl;

    return 0;
}
