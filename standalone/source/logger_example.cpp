#include <iostream>
#include <multiplierless/logger.hpp>

/**
 * @brief Example program demonstrating spdlogger integration
 *
 * This program shows how to use the multiplierless::log_with_spdlog() function
 * to log messages to a file.
 */
int main() {
    std::cout << "Multiplierless Spdlogger Example" << '\n';
    std::cout << "=================================" << '\n';

    // Example 1: Basic logging
    std::cout << "\n[Example 1] Basic logging" << '\n';
    multiplierless::log_with_spdlog("Application started");
    std::cout << "  -> Logged: Application started" << '\n';

    // Example 2: Logging operational messages
    std::cout << "\n[Example 2] Logging operational messages" << '\n';
    multiplierless::log_with_spdlog("Processing data...");
    std::cout << "  -> Logged: Processing data..." << '\n';

    multiplierless::log_with_spdlog("Data processing completed successfully");
    std::cout << "  -> Logged: Data processing completed successfully" << '\n';

    // Example 3: Logging status updates
    std::cout << "\n[Example 3] Logging status updates" << '\n';
    multiplierless::log_with_spdlog("System initialized");
    std::cout << "  -> Logged: System initialized" << '\n';

    multiplierless::log_with_spdlog("Ready to accept connections");
    std::cout << "  -> Logged: Ready to accept connections" << '\n';

    std::cout << "\n=================================" << '\n';
    std::cout << "Check 'multiplierless.log' for logged messages" << '\n';

    return 0;
}
