# AGENTS.md - Agentic Coding Guidelines

This file provides context for AI agents operating in this repository.

## Project Overview

- **Multiplierless FIR Filter optimization** using Ellipsoid Method
- A C++20 library for multiplierless FIR filter design via spectral factorization and convex optimization
- Dependencies: ellalgo-cpp, csd-cpp, ginger-cpp, FFTW, nlohmann_json, spdlog, fmt, doctest, RapidCheck, cxxopts
- Standard: C++20; CMake 3.14+ or xmake

## Build Commands

### xmake (local development, uses sibling checkouts)
```bash
xmake f -m release
xmake -y -j8
xmake test
```
`xmake.lua` expects `../ellalgo-cpp`, `../csd-cpp`, `../ginger-cpp`, and `../lds-cpp`.

### CMake (CI / standalone install)
```bash
cmake -S. -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
```

### Run All Tests
```bash
xmake test
# or via CMake/ctest:
CTEST_OUTPUT_ON_FAILURE=1 ctest --test-dir build --output-on-failure
```

### Run a Single Test
```bash
# xmake
xmake run MultiplierlessTests --test-case="Lowpass Filter*"
# or the built binary / ctest
./build/MultiplierlessTests --test-case="Lowpass Filter*"
ctest --test-dir build -R "MultiplierlessAllTests" -V
```

### Verify Generated Verilog
```bash
python tools/test_verilog.py   # requires iverilog + vvp
```

### Code Formatting (clang-format + cmake-format)
```bash
# C++ sources
clang-format -i include/**/*.hpp source/*.cpp standalone/source/*.cpp test/source/*.cpp
# CMake files (if cmake-format is installed)
cmake-format -i CMakeLists.txt
# Or via CMake's Format.cmake targets
cmake -S. -B build
cmake --build build --target format
cmake --build build --target fix-format
```

### Static Analyzers
```bash
# clang-tidy over the public headers (uses the checks in .clang-tidy)
cmake -S . -B build -DMULTIPLIERLESS_ENABLE_CLANG_TIDY=ON
cmake --build build --target clang-tidy
```

### Sanitizers (optional)
```bash
# Address + Undefined sanitizers
cmake -S. -B build -DUSE_SANITIZER="Address;Undefined"
cmake --build build
```

## Directory Structure

```
multiplierless-cpp/
├── include/multiplierless/    # Public headers (*.hpp)
├── source/                  # Implementation (*.cpp)
├── test/source/             # Test files (*.cpp)
├── standalone/            # Example executable
├── documentation/         # Doxygen config
└── cmake/               # CMake modules
```

## Code Style Guidelines

### Formatting
- **clang-format** with Google-based style (see `.clang-format`)
- Column limit: 100
- Indent width: 4
- Use `BreakBeforeBraces: Attach`

### Naming Conventions
- **Classes**: `CamelCase` (e.g., `LowpassOracle`, `filter_design_construct`)
- **Member variables**: Leading underscore + snake_case (e.g., `_i_Anr`, `_Fdc`)
- **Functions**: `CamelCase` for methods, `snake_case` for free functions
- **Constants**: `kCamelCase` or `UPPER_SNAKE_CASE`

### File Organization
- Header guards: `#pragma once`
- Includes: Sorted (`clang-format -i` will handle)
- Order: related std libs, then external libs, then local headers
- Use explicit namespace qualification (no `using namespace std;`)

### Types
- Use xtensor types: `xt::xarray<double>`, `xt::xtensor<double, N>`
- Use type aliases in classes:
  ```cpp
  using Arr = xt::xarray<double>;
  using Vec = std::valarray<double>;
  ```

### Error Handling
- Return tuples: `std::tuple<Result, bool>` or `std::pair<...>`
- No exceptions in library code unless documented
- Use doctest `CHECK()` / `REQUIRE()` for test assertions

### C++20 Features
- Allowed: designated initializers, structured bindings, `std::optional`, `std::variant`, `std::erase`, `std::set::contains`, `std::string::starts_with`
- Avoid: heap allocations in hot paths
- Use `constexpr` where appropriate

## Testing Guidelines

### Test Framework
- **doctest** for unit tests
- **RapidCheck** for property-based testing (when available)

### Writing Tests
```cpp
#include <doctest/doctest.h>

TEST_CASE("Description") {
    CHECK(condition);
    REQUIRE(condition_that_must_pass);
}
```

### Test Naming
- Use descriptive names: `TEST_CASE("Lowpass Filter (w/ parallel cut)")`
- Group related tests with prefixes

## Common Issues / Gotchas

### Windows MSVC
- Uses `/permissive-` (standards conformance)
- Suppresses C4244 (potential data loss) warnings
- May need `#include <cmath>` for `std::pow`, etc.

### macOS
- Defines `XTENSOR_DISABLE_SVECTOR=1` to avoid Clang template ambiguity with 64-bit `long`/`unsigned long`

### xtensor Usage
- Views are cheap but copies return new arrays
- Use `xt::view()` for slicing without copies
- Use `xt::sum()` with shape inference `()` for scalar result

## Development Workflow

1. Create feature branch: `git checkout -b feature/xxx`
2. Make changes following style guidelines
3. Run formatting: `cmake --build build --target fix-format`
4. Build and test: `cmake --build build && ./build/MultiplierlessTests`
5. Commit: `git add` + `git commit -m "description"`
6. Push and PR

## External Dependencies

CMake uses CPM.cmake (pinned tags) and fetches everything automatically. xmake
uses the sibling checkouts listed above and `add_requires` for FFTW, fmt,
spdlog, nlohmann_json, cxxopts, and doctest.

Key packages:
- `ellalgo-cpp` - ellipsoid method and `Arr`
- `csd-cpp` - CSD conversion
- `ginger-cpp` - root finding
- `FFTW` - FFT-based spectral factorization
- `nlohmann_json` - CLI JSON
- `spdlog` / `fmt` - logging and formatting
- `doctest` / `rapidcheck` - testing