#include <doctest/doctest.h>

#include <multiplierless/verilog_utils.hpp>
#include <stdexcept>
#include <string>
#include <vector>

TEST_CASE("csd_coeffs_to_fixed aligns equal widths") {
    auto const fixed = csd_coeffs_to_fixed({"1.5", "0.25"});
    CHECK(fixed.digits.size() == 2);
    CHECK(fixed.digits[0] == "150");
    CHECK(fixed.digits[1] == "025");
    CHECK(fixed.max_power == 2);
    CHECK(fixed.scale_power == 2);
}

TEST_CASE("csd_coeffs_to_fixed pads integer digits") {
    auto const fixed = csd_coeffs_to_fixed({"0.5", "10.25"});
    CHECK(fixed.digits[0] == "0050");
    CHECK(fixed.digits[1] == "1025");
    CHECK(fixed.max_power == 3);
    CHECK(fixed.scale_power == 2);
}

TEST_CASE("csd_coeffs_to_fixed accepts strings without a decimal point") {
    auto const fixed = csd_coeffs_to_fixed({"101", "1"});
    CHECK(fixed.digits[0] == "101");
    CHECK(fixed.digits[1] == "001");
    CHECK(fixed.max_power == 2);
    CHECK(fixed.scale_power == 0);
}

TEST_CASE("csd_coeffs_to_fixed accepts an empty fraction") {
    auto const fixed = csd_coeffs_to_fixed({"1."});
    CHECK(fixed.digits[0] == "1");
    CHECK(fixed.max_power == 0);
    CHECK(fixed.scale_power == 0);
}

TEST_CASE("csd_coeffs_to_fixed keeps signed digits") {
    auto const fixed = csd_coeffs_to_fixed({"+0-", "-+0"});
    CHECK(fixed.digits[0] == "+0-");
    CHECK(fixed.digits[1] == "-+0");
    CHECK(fixed.max_power == 2);
}

TEST_CASE("csd_coeffs_to_fixed rejects malformed input") {
    CHECK_THROWS_AS(csd_coeffs_to_fixed({}), std::invalid_argument);
    CHECK_THROWS_AS(csd_coeffs_to_fixed({".5"}), std::invalid_argument);
}

TEST_CASE("fix_direct_verilog drops the trailing port comma and rewrites outputs") {
    auto const verilog =
        "module demo (\n"
        "    input clk,\n"
        "    input [7:0] a,\n"
        "    output wire signed [7:0] h0,\n"
        "    output wire signed [7:0] h1,\n"
        ");\n"
        "    wire signed [7:0] h0 = a + 1;\n"
        "    wire signed [7:0] h1 = a - 1;\n"
        "    wire signed [7:0] tmp = a * 2;\n"
        "endmodule\n";

    auto const fixed = fix_direct_verilog(verilog, {"h0", "h1"});
    CHECK(fixed.find("    output wire signed [7:0] h1,\n") == std::string::npos);
    CHECK(fixed.find("    assign h0 = a + 1;\n") != std::string::npos);
    CHECK(fixed.find("    assign h1 = a - 1;\n") != std::string::npos);
    // non-output wires are left alone
    CHECK(fixed.find("    wire signed [7:0] tmp = a * 2;\n") != std::string::npos);
}

TEST_CASE("fix_direct_verilog inserts missing port commas") {
    auto const verilog =
        "module demo (\n"
        "    input clk\n"
        "    output h0\n"
        ");\n"
        "endmodule\n";

    auto const fixed = fix_direct_verilog(verilog, {});
    CHECK(fixed.find("    input clk,\n") != std::string::npos);
    CHECK(fixed.find("    output h0\n") != std::string::npos);
}

TEST_CASE("fix_direct_verilog inserts the comma before a port comment") {
    auto const verilog =
        "module demo (\n"
        "    input clk  // clock\n"
        "    output h0\n"
        ");\n"
        "endmodule\n";

    auto const fixed = fix_direct_verilog(verilog, {});
    CHECK(fixed.find("    input clk  ,// clock\n") != std::string::npos);
}

TEST_CASE("fix_direct_verilog keeps the trailing comma when it is a comment") {
    auto const verilog =
        "module demo (\n"
        "    input clk,  // clock\n"
        "    output h0,\n"
        ");\n"
        "endmodule\n";

    auto const fixed = fix_direct_verilog(verilog, {});
    CHECK(fixed.find("    output h0\n") != std::string::npos);
}

TEST_CASE("fix_direct_verilog passes unrelated text through") {
    CHECK(fix_direct_verilog("", {}) == "");
    CHECK(fix_direct_verilog("module nothing\n", {}) == "module nothing\n");
    // handles CRLF input without treating ');\r' as a port line
    auto const crlf = fix_direct_verilog("module demo (\r\n    input clk\r\n);\r\n", {});
    CHECK(crlf.find("input clk") != std::string::npos);
}

TEST_CASE("fix_direct_verilog only rewrites named output wires") {
    auto const verilog =
        "module demo (\n"
        "    output h0\n"
        ");\n"
        "  wire signed [3:0] h0 = a;\n"
        "  wire signed [3:0] other = b;\n"
        "endmodule\n";

    auto const fixed = fix_direct_verilog(verilog, {"h0"});
    CHECK(fixed.find("  assign h0 = a;\n") != std::string::npos);
    CHECK(fixed.find("  wire signed [3:0] other = b;\n") != std::string::npos);
}
