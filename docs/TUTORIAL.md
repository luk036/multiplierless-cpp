# 🎛️ Multiplierless FIR Filter Design — Tutorial

## From Spec to Silicon: JSON → CSD → Verilog

This tutorial walks you through designing a multiplierless FIR lowpass filter
using the `FirDesign` toolchain. You'll go from a JSON specification file all
the way to a synthesizable Verilog module — ready for ASIC or FPGA.

---

## 📋 Table of Contents

1. [What This Tool Does](#1-what-this-tool-does)
2. [Installation](#2-installation)
3. [Writing Your Filter Spec](#3-writing-your-filter-spec)
4. [Running the CLI](#4-running-the-cli)
5. [Understanding the Output](#5-understanding-the-output)
6. [The Verilog Module](#6-the-verilog-module)
7. [Synthesizing with Yosys](#7-synthesizing-with-yosys)
8. [Full Example Walkthrough](#8-full-example-walkthrough)
9. [Verifying the Filter Meets Specifications](#9-verifying-the-filter-meets-specifications)
10. [Tuning Guide](#10-tuning-guide)

---

## 1. What This Tool Does

```
   JSON spec  →  Ellipsoid optimization  →  CSD quantization  →  Verilog RTL
```

- **Input**: A JSON file describing your desired filter (passband, stopband,
  ripple, attenuation, number of taps, CSD precision).
- **Engine**: The ellipsoid method solves the convex optimization problem
  `min max R(ω)` subject to passband/stopband constraints.
- **Quantization**: Each coefficient is converted to Canonical Signed Digit
  (CSD) format — using only `+1`, `-1`, `0` digits — so hardware
  multiplication becomes shift-and-add.
- **Output**: A JSON file with CSD coefficients AND a synthesizable Verilog
  module with automatic cross-CSE (common sub-expression elimination).

### Why multiplierless?

| Approach | Hardware per tap | 32-tap cost |
|----------|-----------------|-------------|
| Generic multiplier | 1 multiplier | 32 multipliers |
| CSD (nnz=4) | 3 adders/subtractors | 96 adders |
| CSD + Cross-CSE | 2–3 adders (shared) | ~60 adders |

---

## 2. Installation

### Option A: C++ (xmake)

```bash
cd multiplierless-cpp
xmake -y -j8
```

The CLI binary is at `build/windows/x64/release/FirDesign.exe` (or
`build/linux/x86_64/release/FirDesign` on Linux).

### Option B: Python (pip)

```bash
cd ../py/multiplierless
pip install -e .
```

The CLI command is `fir-design`.

### Option C: Python (module)

```bash
conda activate pyflow  # or your environment
python -m multiplierless.fir_design spec.json
```

---

## 3. Writing Your Filter Spec

Create a JSON file, e.g. `my_filter.json`:

```json
{
  "filter_order": 48,
  "passband_edge": 0.12,
  "stopband_edge": 0.20,
  "passband_ripple": 0.025,
  "stopband_attenuation": 0.125,
  "csd_nnz": 6,
  "discretization_factor": 15,
  "max_iters": 50000,
  "tolerance": 1e-14,
  "ellipsoid_radius": 40.0,
  "parallel_cut": true,
  "verilog": {
    "input_width": 16,
    "module_name": "my_fir_filter"
  }
}
```

### Required keys

| Key | Type | Description |
|-----|------|-------------|
| `filter_order` | int (≥4) | Number of FIR taps |
| `csd_nnz` | int (1–16) | Max non-zero CSD digits per coefficient |

All other keys have sensible defaults (shown above).

### Key parameters explained

| Key | Default | What it controls |
|-----|---------|-----------------|
| `passband_edge` | `0.12` | ×π rad. Lower = narrower passband |
| `stopband_edge` | `0.20` | ×π rad. Higher = wider transition band |
| `passband_ripple` | `0.125` | Linear ripple. `0.125` ≈ ±1.0 dB |
| `stopband_attenuation` | `0.125` | Linear attn. `0.125` ≈ −18 dB |
| `discretization_factor` | `15` | Freq samples = N × this. Higher = more accurate but slower |
| `csd_nnz` | *(required)* | Hardware cost = nnz − 1 adders per tap |
| `max_iters` | `50000` | Ellipsoid iteration limit |
| `tolerance` | `1e-14` | Convergence tolerance. Smaller = tighter |
| `ellipsoid_radius` | `40.0` | Initial search region size |
| `parallel_cut` | `true` | Enables faster convergence |
| `spectral_method` | `"fft"` | `"fft"` (Kolmogorov, default, more stable) or `"root"` (Aberth, tunable) |
| `root_tolerance` | `1e-8` | Aberth convergence tolerance (only used when `spectral_method="root"`). Larger = looser but faster |

### Spectral Factorization Methods

The CLI supports two spectral factorization algorithms:

| Method | Key | Speed | Accuracy | When to use |
|--------|-----|-------|----------|-------------|
| **FFT (Kolmogorov 1939)** | `"fft"` | Stable, slower (~1482 iter) | Reference (exact values) | **Default**, most stable |
| **Aberth root-finding** | `"root"` | Fast (~333 iter) | Good (configurable tolerance) | Quick experiments, tuning |

The FFT method is the **default** because it's numerically more stable —
especially for large filters and tight specifications. The root-finding method
can be faster but requires tuning `root_tolerance` for each filter order:
- `1e-4`: Loose but fast — good for quick experiments
- `1e-8`: Balanced speed and accuracy
- `1e-12`: Tight — use for narrow transition bands or high-order filters

> **Note**: For very large filters (N > 64), use the default `"fft"` method.
> Switch to `"root"` only when you need faster iteration counts and are
> willing to tune `root_tolerance`.

### To get Verilog output

Add a `verilog` section:

```json
"verilog": {
    "input_width": 16,
    "module_name": "fir_filter",
    "form": "transpose"
  }
}```

Omit the `verilog` key entirely to get coefficients only (faster).

### Verilog form: direct vs transpose

| Form | Description | Output |
|------|-------------|--------|
| `"transpose"` **(default)** | Complete pipelined FIR filter with clock/reset | Single `y` = sum(hᵢ · x) |
| `"direct"` | Multiplier bank — each tap computed in parallel | N separate `hᵢ` = coeffᵢ · x |

The transpose form is a **complete FIR filter module** with built-in pipeline
registers — instantiate and connect. The direct form provides raw coefficient
multiplier outputs for custom integration.

> **Note**: The `form` option only matters if you include a `verilog` section.
> Omit the `verilog` key entirely for coefficient-only output.

---

## 4. Running the CLI

### C++

```bash
./build/windows/x64/release/FirDesign.exe my_filter.json | tee output.json
```

### Python

```bash
fir-design my_filter.json | tee output.json
```

The JSON output is printed to stdout. Redirect to save:

```bash
fir-design my_filter.json > output.json
```

**Important**: The optimization can take 1–5 minutes depending on `filter_order`
and `max_iters`. For first experiments, try `{"filter_order": 16, "csd_nnz": 4}`
which converges in seconds.

---

## 5. Understanding the Output

```json
{
  "filter_order": 48,
  "csd_nnz": 6,
  "iterations": 2341,
  "coefficients": [
    {"index": 0,  "value": 0.008912, "csd": "0.000000+00+000-000-0+0-"},
    {"index": 1,  "value": 0.017234, "csd": "0.00000+000+0+00-0+00-"},
    ...
    {"index": 47, "value": 0.002145, "csd": "0.00000000+000-0+0+"}
  ],
  "verilog": "module fir_filter( ..."
}
```

### Output fields

| Field | Meaning |
|-------|---------|
| `iterations` | Ellipsoid method iterations performed. Fewer = faster convergence |
| `coefficients[].value` | CSD-quantized coefficient (double). Use in simulation |
| `coefficients[].csd` | CSD string. `+` = +1, `-` = −1, `0` = 0, `.` = decimal point |
| `verilog` | Full synthesizable Verilog module (if requested) |

### Reading the CSD string

```
0.000+00-0+0+
│ │   │  │ │ │
│ │   +  - + +  ← non-zero digits (4 used out of nnz=6)
│ └────────────  fractional part (powers: 2⁻¹, 2⁻², ...)
└──────────────  integer part (powers: 2⁰, 2¹, 2², ...)
```

Each non-zero digit represents one adder/subtractor in hardware.
The total adder count per tap = `nnz − 1` at most.

### Filter symmetry

The impulse response `h` is **minimum-phase** (from spectral factorization).
For a linear-phase FIR filter, use the squared magnitude response
(or autocorrelation) directly.

---

## 6. The Verilog Module

The generated Verilog uses **only shift and add/subtract** operations — no `*`.

By default the generator produces a **transpose-form** FIR filter with built-in
pipeline registers. You can select the **direct-form** multiplier bank via
`"form": "direct"` in the JSON spec.

### Transpose-form (default)

The transpose-form module is a **complete FIR filter** ready for instantiation:

### Module signature

```verilog
module fir_filter (
    input clk,                     // Clock
    input rst_n,                   // Active-low reset
    input signed [15:0] x,         // 16-bit input sample
    output signed [41:0] y         // Filter output = Σ hᵢ · x[n−i]
);
```

### Internal structure

```verilog
// Shifted copies of input (one wire per unique shift amount)
wire signed [41:0] x_shift17 = x <<< 17;
wire signed [41:0] x_shift16 = x <<< 16;
wire signed [41:0] x_shift14 = x <<< 14;
...

// Cross-CSE: shared sub-expression
// Pattern "+0-" found in multiple coefficients → factored out!
wire signed [41:0] _cse_0 = x_shift14 - x_shift12;

// Transpose-form pipeline registers
reg signed [41:0] sum0;
reg signed [41:0] sum1;
...
reg signed [41:0] sum47;

// Coefficients applied in reverse order (canonical transpose form)
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        sum0 <= 0; sum1 <= 0; ... sum47 <= 0;
    end else begin
        sum0 <= /* h[47]*x using CSD */;               // last coefficient
        sum1 <= sum0 + /* h[46]*x using CSD */;
        ...
        sum47 <= sum46 + /* h[0]*x using CSD */;        // first coefficient
    end
end

assign y = sum47;
```

### Direct-form (opt-in)

Set `"form": "direct"` in the JSON spec to generate a **multiplier bank**
with separate tap outputs (for custom external accumulation):

### Module signature (direct-form)

```verilog
module fir_filter (
    input signed [15:0] x,       // 16-bit input sample
    output signed [41:0] h0,     // coefficient h[0] (41 bits = 16 + 25)
    output signed [41:0] h1,
    ...
    output signed [41:0] h47
);
```

### Internal structure (direct-form)

```verilog
// Shifted copies of input (one wire per unique shift amount)
wire signed [41:0] x_shift17 = x <<< 17;
wire signed [41:0] x_shift16 = x <<< 16;
wire signed [41:0] x_shift14 = x <<< 14;
...

// Cross-CSE: shared sub-expression
wire signed [41:0] _cse_0 = x_shift14 - x_shift12;

// Per-coefficient assignment (shift-and-add only)
wire signed [41:0] h0 = (_cse_0 >>> 6) + x_shift4 + x_shift0;
wire signed [41:0] h5 = (_cse_0 >>> 5) + x_shift3 + x_shift0;
...
```

### Cross-CSE explained

When multiple coefficients share the same CSD pattern (e.g., `+0-`), it is
computed once as a shared wire `_cse_0` and reused across taps with
right-shifts (`>>>`). This **reduces total adder count** across the entire
filter.

```
Without CSE:  h0 needs 1 adder for +0-, h5 needs 1 adder for +0- = 2 adders
With CSE:     _cse_0 = +0- (1 adder), h0 = _cse_0 >>> 6, h5 = _cse_0 >>> 5 = 1 adder
```

### Output width

Output width = `input_width + max_power`, where `max_power` is the highest
shift amount across all coefficients.

| input_width | Typical output_width |
|-------------|---------------------|
| 8 | ~25 |
| 16 | ~41 |
| 24 | ~55 |

---

## 7. Synthesizing with Yosys

The Verilog output is ready for Yosys synthesis:

```bash
# Synthesize to a generic gate-level netlist
yosys -p "
  read_verilog fir_filter_cpp.v;
  synth -top fir_filter;
  write_json fir_filter_synth.json;
  stat
"
```

To target a specific technology (e.g., SkyWater 130nm):

```bash
yosys -p "
  read_verilog fir_filter_cpp.v;
  synth -top fir_filter;
  abc -liberty /path/to/sky130_fd_sc_hd__tt_025C_1v80.lib;
  write_verilog fir_filter_mapped.v;
  stat
"
```

### Expected area

For a 32-tap filter with nnz=4 and cross-CSE, the actual Yosys synthesis
results from our example are:

```
=== fir_filter ===

   Number of wires:                100
   Number of wire bits:           3778
   Number of public wires:          49
   Number of cells:                 82
     $add                           51
     $neg                            7
     $sub                           24

   No multipliers — purely shift-and-add!
```

Without cross-CSE: ~110 cells. With cross-CSE: **82 cells** (a 25% reduction
from shared sub-expression elimination). Each `$add` or `$sub` is one
hardware adder/subtractor.

To run synthesis with your own filter, use the provided script:

```bash
# Fix Verilog port syntax (CSD generator produces comma-less ports)
python tools/fix_verilog_ports.py fir_filter_tutorial.v

# Synthesize with Yosys
yosys -Q synthesize.ys
```

The `synthesize.ys` script does the following:
```tcl
read_verilog fir_filter_tutorial.v    # load your design
hierarchy -top fir_filter             # check module hierarchy
proc                                  # process always blocks
opt                                   # optimize expressions
stat                                  # print area/adder counts
write_json fir_filter_synth.json      # JSON netlist for analysis
write_verilog -noattr fir_filter_synth.v  # synthesized Verilog
```

---

## 8. Full Example Walkthrough

Let's design a 32-tap lowpass filter with 4 non-zero CSD digits per
coefficient.

### Step 1: Write the spec

Save as `lowpass_32.json`:

```json
{
  "filter_order": 32,
  "csd_nnz": 4,
  "passband_edge": 0.12,
  "stopband_edge": 0.20,
  "passband_ripple": 0.125,
  "stopband_attenuation": 0.125,
  "verilog": {
    "input_width": 16,
    "module_name": "fir_filter"
  }
}
```

### Step 2: Run the tool

```bash
fir-design lowpass_32.json | tee lowpass_32_output.json
```

### Step 3: Check the results

Look at the output:

```json
{
  "iterations": 1482,
  "filter_order": 32,
  "csd_nnz": 4,
  "coefficients": [
    {"index": 0,  "value": 0.008724, "csd": "0.000000+00+000-000-"},
    {"index": 7,  "value": 0.128059, "csd": "0.00+0000+0-00+"},
    {"index": 15, "value": -0.017375,"csd": "0.00000-00-00+0-"},
    ...
  ]
}
```

- **1482 iterations** — converged in < 2 minutes
- **Coefficient h[7]** = 0.128059, CSD = `0.00+0000+0-00+` — 4 non-zero digits
- **Coefficient h[15]** = −0.017375, CSD = `0.00000-00-00+0-` — 4 non-zero digits

### Step 4: Extract the Verilog

The `verilog` field contains the full module. To extract it to a file:

```bash
python -c "
import json, sys
data = json.load(open('lowpass_32_output.json'))
with open('fir_filter.v', 'w') as f:
    f.write(data['verilog'])
print(f'{len(data[\"verilog\"])} chars written')
"
```

Or use the included tool:

```bash
python tools/extract_verilog.py
```

### Step 5: Fix Verilog syntax & synthesize

The CSD generator outputs ports without commas between them. Fix with:

```bash
python tools/fix_verilog_ports.py fir_filter_tutorial.v
```

Then synthesize:

```bash
yosys -Q synthesize.ys
```

Expected output:
```
=== fir_filter ===
   Number of cells:   82
     $add             51
     $neg              7
     $sub             24
```

### Step 6: Integrate into your design

Instantiate the module in your top-level:

```verilog
wire signed [15:0] sample_in;
wire signed [41:0] tap_0, tap_1, ..., tap_31;

fir_filter u_fir (
    .x(sample_in),
    .h0(tap_0), .h1(tap_1), ..., .h31(tap_31)
);

// Sum the taps with appropriate delays to get the filter output
// (or feed into a shift-register + adder tree)
```

---

## 9. Verifying the Filter Meets Specifications

After generating coefficients, verify the filter actually meets your
passband/stopband requirements before committing to silicon.

### 9.1 Frequency Response Analysis

Use the output JSON to compute the frequency response in Python:

```python
import json
import numpy as np

# Load the design output
with open("lowpass_32_output.json") as f:
    data = json.load(f)

h = np.array([c["value"] for c in data["coefficients"]])
N = len(h)

# Compute frequency response
w = np.linspace(0, np.pi, 4096)  # fine grid
H = np.zeros(len(w), dtype=complex)
for i, wi in enumerate(w):
    H[i] = sum(h[k] * np.exp(-1j * k * wi) for k in range(N))

# Magnitude squared response (what the optimizer works with)
R = np.abs(H) ** 2

# Find indices for passband and stopband edges
wpass = int(0.12 * len(w))   # passband: 0 to 0.12π
wstop = int(0.20 * len(w))   # stopband: 0.20π to π

# Passband check: should be within [L², U²]
# (values from the JSON spec: ripple=0.125 → L² ≈ 0.64, U² ≈ 1.44)
passband_ok = np.all(R[:wpass] >= 0.64) and np.all(R[:wpass] <= 1.44)

# Stopband check: should be below attenuation target
# (from spec: attenuation=0.125 → S² ≈ 0.0156)
stopband_max = np.max(R[wstop:])
stopband_ok = stopband_max <= 0.0156

print(f"Passband OK: {passband_ok}")
print(f"Stopband max R(ω): {stopband_max:.6f}")
print(f"Stopband OK: {stopband_ok}")
```

### 9.2 Verify with scipy.signal.freqz

A simpler approach using SciPy (if installed):

```python
from scipy.signal import freqz
import numpy as np
import json
import matplotlib.pyplot as plt

with open("lowpass_32_output.json") as f:
    data = json.load(f)

h = np.array([c["value"] for c in data["coefficients"]])
N = len(h)

w, H = freqz(h, worN=4096)
w = w / np.pi  # normalize to [0, 1]

# Plot
plt.figure(figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.plot(h, 'o-', markersize=3)
plt.title(f"Impulse Response ({N} taps)")
plt.xlabel("Tap index")
plt.ylabel("h[n]")

plt.subplot(1, 2, 2)
plt.plot(w, 20 * np.log10(np.abs(H) + 1e-12))
plt.axvline(0.12, color='orange', linestyle='--', label='passband edge')
plt.axvline(0.20, color='red', linestyle='--', label='stopband edge')
plt.ylim(-80, 5)
plt.title("Magnitude Response (dB)")
plt.xlabel("Normalized frequency (×π)")
plt.ylabel("|H(ω)| [dB]")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("fir_response.png", dpi=150)
plt.show()
```

### 9.3 What to check

| Check | Expected | How |
|-------|----------|-----|
| **Passband flatness** | R(ω) within [L², U²] for ω ∈ [0, ωₚ] | `np.min(R[:wpass])` ≥ L², `np.max(R[:wpass])` ≤ U² |
| **Stopband attenuation** | R(ω) ≤ S² for ω ∈ [ωₛ, π] | `np.max(R[wstop:])` ≤ S² |
| **Impulse response shape** | Symmetric about center (linear phase) or minimum-phase | Visual check of h coefficients |
| **No NaN/Inf** | All coefficients finite | `np.all(np.isfinite(h))` |
| **CSD digit count** | ≤ csd_nnz per coefficient | Count `+` and `-` in CSD string |

### 9.4 Verification script

The project includes a verification helper:

```bash
python tools/verify_filter.py lowpass_32_output.json
```

Output shows pass/fail for each specification band.

### 9.5 If the filter fails verification

The verification script checks at **two resolutions**:

```
--- Oracle grid (480 points, disc_factor=15) ---
[PASS] Passband: R(w) in [0.7932, 1.2646]
[PASS] Stopband: max R(w) = 0.000000 <= 0.015625

--- Fine grid (16384 points, advisory) ---
[FAIL] Passband: R(w) in [0.7679, 1.2647]
       (discretization artifact: increase disc_factor)
[PASS] Stopband: max R(w) = 0.000000 <= 0.015625

ALL CHECKS PASSED
```

**The oracle grid is authoritative**: the optimizer guarantees feasibility
only at its sampled frequency points (15×N). The fine-grid passband dip
(0.7679 vs 0.7901) occurs **between** the oracle's sampling points — it's a
discretization artifact, not a quantization error.

> The discrete optimizer (`cutting_plane_optim_q`) with CSD quantization returns
> `rcsd` — the quantized autocorrelation — which meets specs at every sampled
> frequency. No quantization error is present.

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| Fail at oracle grid (rare) | csd_nnz too low or bad convergence | Increase `max_iters` or `csd_nnz` |
| Fail at fine grid only | Discretization gap between sampled points | Increase `discretization_factor` to 30 |
| Stopband attenuation too weak | Not enough taps | Increase `filter_order` by 8–16 |
| Erratic frequency response | Convergence tolerance too loose | Decrease `tolerance` to `1e-20` |

---

## 10. Tuning Guide

### "My filter doesn't converge"

Try in order:
1. **Decrease `passband_ripple`** — larger ripple (e.g. 0.2) is easier to achieve
2. **Increase `max_iters`** — 50000 → 100000 gives the solver more time
3. **Loosen `tolerance`** — `1e-14` → `1e-10` stops earlier
4. **Decrease `filter_order`** — fewer taps = smaller problem

### "I want fewer adders (smaller area)"

- **Decrease `csd_nnz`** — from 4 → 3 saves ~25% adders but may degrade
  stopband performance. Check the output coefficients to see if stopband
  attenuation is still acceptable.

### "I need better stopband attenuation"

- **Increase `filter_order`** — more taps = sharper roll-off
- **Increase `csd_nnz`** — more non-zero digits = better approximation of
  ideal coefficients
- **Narrow the transition band** — increase `stopband_edge` closer to
  `passband_edge` (but this makes convergence harder)

### "I need more accurate frequency response"

- **Increase `discretization_factor`** — from 15 → 30 doubles the frequency
  sampling points (at the cost of more memory and slower iterations)

### Quick reference: parameter impact

| Parameter | Increase → | Effect |
|-----------|-----------|--------|
| `filter_order` | More taps | Better stopband, slower |
| `csd_nnz` | More digits | Better accuracy, more area |
| `max_iters` | More iterations | Better convergence, slower |
| `tolerance` | Larger value | Faster exit, looser optimality |
| `discretization_factor` | More samples | More accurate constraints, more memory |
| `passband_ripple` | Larger | Easier to design, more ripple |
| `stopband_attenuation` | Smaller | Deeper stopband, harder to converge |
| `spectral_method` | `"root"` | Use root-finding, faster iterations but needs tuning |
| `root_tolerance` | Larger value | Faster convergence, looser root accuracy |

---

## 🔧 Troubleshooting

| Problem | Solution |
|---------|----------|
| `Error: cannot open file` | Check path — use absolute or relative to CWD |
| `JSON parse error` | Validate JSON at [jsonlint.com](https://jsonlint.com) |
| `Optimization failed — no feasible solution` | See tuning guide above |
| CLI hangs | Reduce `max_iters` or `filter_order` for testing |
| `fir-design: command not found` | Run `pip install -e .` in the Python project |
| No Verilog in output | Add `"verilog": {"input_width": 16}` to your JSON spec |
| "Aberth did not converge" | Increase `root_tolerance` (e.g. `1e-6`) or switch to `"fft"` |

---

## 📚 Reference

- **Paper**: Wu, Boyd & Vandenberghe — "FIR Filter Design via Spectral
  Factorization and Convex Optimization"
- **C++ repo**: [luk036/multiplierless-cpp](https://github.com/luk036/multiplierless-cpp)
- **Python repo**: [luk036/multiplierless](https://github.com/luk036/multiplierless)
- **CSD library (C++)**: [luk036/csd-cpp](https://github.com/luk036/csd-cpp)
- **CSD library (Python)**: [luk036/csdigit](https://github.com/luk036/csdigit)
- **Ellipsoid method**: [luk036/ellalgo](https://github.com/luk036/ellalgo)
