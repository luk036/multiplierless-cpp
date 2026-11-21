[![Actions Status](https://github.com/luk036/multiplierless-cpp/workflows/macOS/badge.svg)](https://github.com/luk036/multiplierless-cpp/actions)
[![Actions Status](https://github.com/luk036/multiplierless-cpp/workflows/Windows/badge.svg)](https://github.com/luk036/multiplierless-cpp/actions)
[![Actions Status](https://github.com/luk036/multiplierless-cpp/workflows/Ubuntu/badge.svg)](https://github.com/luk036/multiplierless-cpp/actions)
[![Actions Status](https://github.com/luk036/multiplierless-cpp/workflows/Install/badge.svg)](https://github.com/luk036/multiplierless-cpp/actions)
[![codecov](https://codecov.io/gh/luk036/multiplierless-cpp/branch/master/graph/badge.svg)](https://codecov.io/gh/luk036/multiplierless-cpp)

# ♒︎ Multiplierless FIR Filter Design (C++)

A C++20 library and CLI for designing multiplierless FIR lowpass filters with
the ellipsoid method. Coefficients are quantized to Canonical Signed Digit
(CSD) form, so each tap is implemented with shifts and a small number of
adders/subtractors instead of a general-purpose multiplier.

The optimizer works on the squared-magnitude response
`R(ω) = r₀ + 2·Σ rₖ·cos(kω)`, then recovers a minimum-phase impulse response by
spectral factorization (FFT or Aberth root-finding). The CLI can emit a
synthesizable Verilog module in transpose or direct form with cross-CSE.

## Build

### xmake (recommended for local development)

`xmake.lua` builds against the sibling checkouts `../ellalgo-cpp`,
`../csd-cpp`, `../ginger-cpp`, and `../lds-cpp`:

```bash
xmake f -m release
xmake -y -j8
xmake test
```

The CLI is `build/windows/x64/release/FirDesign.exe` on Windows and
`build/linux/x86_64/release/FirDesign` on Linux.

### CMake

CMake fetches pinned releases through CPM.cmake and is what the GitHub
workflows use:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

FFTW, OpenBLAS, and the C++ toolchain must be available. On Windows the
workflow installs them with conda.

## Run

```bash
FirDesign sample_filter.json
FirDesign sample_filter.json -o result.json --verilog-out fir_filter.v
FirDesign sample_filter_direct.json -o direct.json --verilog-out fir_filter_direct.v
```

The input schema is `fir_design_input.schema.json`; the output schema is
`fir_design_output.schema.json`.

### Output

Each coefficient is reported as an exact pair:

```json
{"index": 0, "value": 0.01275634765625, "csd": "0.00000+0-0+000+"}
```

`value` is derived from the emitted CSD string, so the numeric and string
views cannot disagree. All taps share one right-aligned fixed-point scale
when handed to the Verilog generator.

### Verilog forms

- `"form": "transpose"` (default): a complete pipelined FIR filter with
  `clk`/`rst_n` and one `y` output. The accumulator includes
  `ceil(log2(taps))` bits of headroom.
- `"form": "direct"`: a multiplier bank with one output per tap (`h0..hN`).

## Verify

```bash
xmake test                    # doctest + RapidCheck suite
python tools/test_verilog.py  # compile + impulse-simulate both Verilog forms
```

`tools/test_verilog.py` requires Icarus Verilog (`iverilog`/`vvp`) and checks
the coefficient/CSD consistency, fixed-point alignment, compilation, and an
impulse-response simulation for both forms.

## Dependencies

- [ellalgo-cpp](https://github.com/luk036/ellalgo-cpp) — ellipsoid method / `Arr`
- [csd-cpp](https://github.com/luk036/csd-cpp) — CSD conversion
- [ginger-cpp](https://github.com/luk036/ginger-cpp) — root finding
- [FFTW](https://www.fftw.org/) — FFT-based spectral factorization
- [nlohmann/json](https://github.com/nlohmann/json) — CLI JSON
- [spdlog](https://github.com/gabime/spdlog) / [fmt](https://github.com/fmtlib/fmt) — logging
- [doctest](https://github.com/doctest/doctest) / [RapidCheck](https://github.com/emil-e/rapidcheck) — tests

## Documentation

- [Tutorial](docs/TUTORIAL.md) — JSON → CSD → Verilog walkthrough
- [中文教程](docs/TUTORIAL_CN.md)
- [Doxygen documentation](https://luk036.github.io/multiplierless-cpp)

## License

MIT. See [LICENSE](LICENSE).
