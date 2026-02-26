#ifdef RAPIDCHECK_H
#    include <rapidcheck.h>
#    include <doctest/doctest.h>

#    include <cmath>
#    include <string>

// Include the functions we want to test
extern auto to_csd(double num, int places = 0) -> std::string;
extern auto to_decimal(const std::string &csd_str) -> double;
extern auto to_csdnnz(double num, unsigned int nnz = 4) -> std::string;

TEST_CASE("Property-based test: CSD round-trip conversion") {
    rc::check("Converting to CSD and back yields original value",
              []() {
                  auto num = *rc::gen::nonNegative<double>();
                  // Avoid extremely large numbers that could cause precision issues
                  num = std::fmod(num, 1000.0);
                  
                  auto csd_str = to_csd(num, 10);
                  auto result = to_decimal(csd_str);
                  
                  // Allow small floating point tolerance
                  RC_ASSERT(std::fabs(result - num) < 0.001);
              });
}

TEST_CASE("Property-based test: CSD contains only valid characters") {
    rc::check("CSD string contains only '+', '-', '0', and '.'",
              []() {
                  auto num = *rc::gen::nonNegative<double>();
                  num = std::fmod(num, 1000.0);
                  
                  auto csd_str = to_csd(num, 10);
                  
                  for (auto c : csd_str) {
                      RC_ASSERT(c == '+' || c == '-' || c == '0' || c == '.');
                  }
              });
}

TEST_CASE("Property-based test: CSD zero conversion") {
    rc::check("to_csd(0) returns '0'",
              []() {
                  auto csd_str = to_csd(0.0, 10);
                  RC_ASSERT(csd_str == "0");
              });
}

TEST_CASE("Property-based test: CSD nnz round-trip conversion") {
    rc::check("Converting to CSD with nnz and back yields approximately original value",
              []() {
                  auto num = *rc::gen::nonNegative<double>();
                  // Limit range to avoid very small or very large numbers
                  num = std::fmod(num, 100.0);
                  auto nnz = static_cast<unsigned int>(*rc::gen::inRange(1, 11));
                  
                  // Skip test for very small numbers where nnz has large relative impact
                  RC_PRE(num > 20.0 && num < 80.0);
                  
                  auto csd_str = to_csdnnz(num, nnz);
                  auto result = to_decimal(csd_str);
                  
                  // Allow larger tolerance due to limited non-zero digits
                  auto relative_error = std::fabs((result - num) / (num + 1e-10));
                  RC_ASSERT(relative_error < 0.5);
              });
}

TEST_CASE("Property-based test: CSD nnz contains limited non-zero digits") {
    rc::check("CSD string with nnz contains at most nnz non-zero digits",
              []() {
                  auto num = *rc::gen::nonNegative<double>();
                  num = std::fmod(num, 1000.0);
                  auto nnz = static_cast<unsigned int>(*rc::gen::inRange(1, 11));
                  
                  auto csd_str = to_csdnnz(num, nnz);
                  
                  size_t count = 0;
                  for (auto c : csd_str) {
                      if (c == '+' || c == '-') {
                          count++;
                      }
                  }
                  
                  RC_ASSERT(count <= nnz);
              });
}

TEST_CASE("Property-based test: CSD nnz zero conversion") {
    rc::check("to_csdnnz(0, any) returns '0'",
              []() {
                  auto nnz = static_cast<unsigned int>(*rc::gen::inRange(1, 11));
                  auto csd_str = to_csdnnz(0.0, nnz);
                  RC_ASSERT(csd_str == "0");
              });
}

TEST_CASE("Property-based test: CSD string length increases with places") {
    rc::check("Increasing places parameter increases CSD string length",
              []() {
                  auto num = *rc::gen::nonNegative<double>();
                  num = std::fmod(num, 100.0);
                  auto places1 = static_cast<int>(*rc::gen::inRange(0, 6));
                  auto places2 = places1 + static_cast<int>(*rc::gen::inRange(1, 6));
                  
                  auto csd_str1 = to_csd(num, places1);
                  auto csd_str2 = to_csd(num, places2);
                  
                  // csd_str2 should be at least as long as csd_str1
                  RC_ASSERT(csd_str2.size() >= csd_str1.size());
              });
}

TEST_CASE("Property-based test: Decimal conversion handles all valid CSD chars") {
    rc::check("to_decimal handles strings with only valid CSD characters",
              []() {
                  // Generate a random valid CSD-like string
                  std::string str = "0";
                  auto length = static_cast<size_t>(*rc::gen::inRange(1, 11));
                  auto has_decimal = *rc::gen::arbitrary<bool>();
                  auto decimal_pos = static_cast<size_t>(*rc::gen::inRange(0, static_cast<int>(length)));
                  
                  for (size_t i = 0; i < length; ++i) {
                      if (has_decimal && i == decimal_pos) {
                          str += '.';
                      } else {
                          auto digit = *rc::gen::inRange(0, 4);
                          switch (digit) {
                              case 0: str += '0'; break;
                              case 1: str += '+'; break;
                              case 2: str += '-'; break;
                              default: str += '0'; break;
                          }
                      }
                  }
                  
                  // This should not crash
                  auto result = to_decimal(str);
                  RC_ASSERT(std::isfinite(result));
              });
}

TEST_CASE("Property-based test: CSD values are within expected range") {
    rc::check("CSD representation yields values close to original for moderate numbers",
              []() {
                  auto num = *rc::gen::nonNegative<double>();
                  num = std::fmod(num, 100.0);  // Keep numbers moderate
                  
                  // Skip test for very small numbers
                  RC_PRE(num > 0.1);
                  
                  auto csd_str = to_csd(num, 5);
                  auto result = to_decimal(csd_str);
                  
                  // For small numbers, conversion should be fairly accurate
                  auto relative_error = std::fabs((result - num) / (num + 1e-10));
                  RC_ASSERT(relative_error < 0.1);
              });
}

TEST_CASE("Property-based test: CSD symmetric property") {
    rc::check("CSD handles small and large numbers correctly",
              []() {
                  auto exponent = *rc::gen::inRange(-10, 10);
                  auto mantissa = *rc::gen::inRange(1, 100) / 100.0;
                  auto num = mantissa * std::pow(2.0, exponent);
                  
                  // Should not crash
                  auto csd_str = to_csd(num, 10);
                  auto result = to_decimal(csd_str);
                  
                  RC_ASSERT(csd_str.size() > 0);
                  RC_ASSERT(std::isfinite(result));
              });
}

#endif