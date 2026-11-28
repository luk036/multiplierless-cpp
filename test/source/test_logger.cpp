#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <multiplierless/logger.hpp>
#include <string>

TEST_CASE("log_with_spdlog writes to the log file") {
    auto const log_file = std::filesystem::path{"multiplierless.log"};
    std::filesystem::remove(log_file);

    multiplierless::log_with_spdlog("unit test message");

    REQUIRE(std::filesystem::exists(log_file));
    std::ifstream in(log_file);
    std::string line;
    bool found = false;
    while (std::getline(in, line)) {
        if (line.find("unit test message") != std::string::npos) found = true;
    }
    CHECK(found);
}
