# AGENTS.md - Agentic Coding Guidelines

This file provides context for AI agents operating in this repository.

## Project Overview

- **Multiplierless FIR Filter optimization** using Ellipsoid Method
- A C++17 library for multiplierless FIR filter design via spectral factorization and convex optimization
- Dependencies: xtensor, spdlog, fmt, doctest, RapidCheck (property-based testing)
- Standard: C++17, CMake 3.14+

## Build Commands

### Full Build
```bash
cmake -S. -B build
cmake --build build
```

### Run All Tests
```bash
./build/test/MultiplierlessTests
# or via ctest:
cd build/test && CTEST_OUTPUT_ON_FAILURE=1 ctest
```

### Run a Single Test
```bash
# Option 1: Use doctest filter with command-line argument
./build/test/MultiplierlessTests --test-case="Lowpass Filter*"
# Option 2: Use ctest with test name
cd build/test && ctest -R "multiplierlessTests" -V
```

### Code Formatting (clang-format + cmake-format)
```bash
cmake -S. -B build
# View changes without applying
cmake --build build --target format
# Apply fixes
cmake --build build --target fix-format
```

### Static Analysis (optional)
```bash
# Enable clang-tidy during build
cmake -S. -B build -DUSE_STATIC_ANALYZER=clang-tidy
cmake --build build
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

### C++17 Features
- Allowed: structured bindings (`auto [a, b] = ...;`), `std::optional`, `std::variant`
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
4. Build and test: `cmake --build build && ./build/test/MultiplierlessTests`
5. Commit: `git add` + `git commit -m "description"`
6. Push and PR

## External Dependencies (CPM.cmake)

Managed via CPM.cmake - fetched automatically on first build. Key packages:
- `xtensor` - multi-dimensional arrays
- `spdlog` - logging
- `fmt` - formatting
- `doctest` - testing
- `rapidcheck` - property-based testing