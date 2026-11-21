"""Run FirDesign and verify the generated Verilog with Icarus Verilog.

Checks, for both transpose and direct forms:
  * every reported ``value`` equals the value encoded by its ``csd`` string;
  * every tap is right-aligned to one shared fixed-point scale;
  * the generated module compiles and passes an impulse-response testbench.

Usage:
    python tools/test_verilog.py                       # both forms, default specs
    python tools/test_verilog.py --cli path/to/FirDesign
"""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

IVERILOG = shutil.which("iverilog")
VVP = shutil.which("vvp")
ROOT = Path(__file__).resolve().parent.parent


def parse_csd(csd: str) -> float:
    int_part, _, frac_part = csd.partition(".")
    value = 0.0
    for i, ch in enumerate(int_part):
        exp = len(int_part) - 1 - i
        if ch == "+":
            value += 2.0**exp
        elif ch == "-":
            value -= 2.0**exp
    for i, ch in enumerate(frac_part):
        exp = -(i + 1)
        if ch == "+":
            value += 2.0**exp
        elif ch == "-":
            value -= 2.0**exp
    return value


def align(csd_list: list[str]) -> tuple[list[str], int]:
    split = []
    for csd in csd_list:
        int_part, _, frac_part = csd.partition(".")
        split.append((int_part, frac_part))
    int_width = max(len(ip) for ip, _ in split)
    frac_width = max(len(fp) for _, fp in split)
    raws = [ip.rjust(int_width, "0") + fp.ljust(frac_width, "0") for ip, fp in split]
    return raws, frac_width


def raw_to_int(raw: str) -> int:
    value = 0
    for i, ch in enumerate(raw):
        power = len(raw) - 1 - i
        if ch == "+":
            value += 1 << power
        elif ch == "-":
            value -= 1 << power
    return value


def output_width(verilog: str, port: str) -> int:
    marker = f"output signed ["
    for line in verilog.splitlines():
        if marker in line and port in line:
            width = line.split("[", 1)[1].split(":", 1)[0]
            return int(width) + 1
    raise ValueError(f"cannot find output width for {port}")


def signed_literal(value: int, width: int) -> str:
    return f"-{width}'sd{abs(value)}" if value < 0 else f"{width}'sd{value}"


def gen_transpose_tb(payload: dict, expected: list[int], width: int) -> str:
    n = len(expected)
    cases = "\n".join(
        f"            {i}: expected = {signed_literal(e, width)};"
        for i, e in enumerate(expected)
    )
    return f"""`timescale 1ns/1ps
module tb;
    reg clk, rst_n;
    reg signed [15:0] x;
    wire signed [{width - 1}:0] y;
    fir_filter uut(.clk(clk), .rst_n(rst_n), .x(x), .y(y));
    initial clk = 0;
    always #5 clk = ~clk;
    integer i, errors;
    reg signed [{width - 1}:0] expected;
    initial begin
        errors = 0; rst_n = 0; x = 0; #15; rst_n = 1; #5;
        for (i = 0; i < {n}; i = i + 1) begin
            x = (i == 0) ? 1 : 0;
            @(posedge clk); #1;
            case (i)
{cases}
            endcase
            if (y !== expected) begin
                $display("FAIL cycle %0d: y=%0d expected=%0d", i, y, expected);
                errors = errors + 1;
            end
        end
        if (errors == 0) $display("ALL TESTS PASSED");
        else $display("FAILED: %0d errors", errors);
        $finish;
    end
endmodule
"""


def gen_direct_tb(payload: dict, expected: list[int], widths: list[int]) -> str:
    n = len(expected)
    ports = "\n".join(f"    wire signed [{widths[i] - 1}:0] h{i};" for i in range(n))
    conns = ",\n".join(f"        .h{i}(h{i})" for i in range(n))
    checks = "\n".join(
        f"        if (h{i} !== {signed_literal(expected[i], widths[i])}) begin\n"
        f'            $display("FAIL h{i}: %0d expected %0d", h{i}, {expected[i]});\n'
        f"            errors = errors + 1;\n"
        f"        end"
        for i in range(n)
    )
    return f"""`timescale 1ns/1ps
module tb;
    reg signed [15:0] x;
{ports}
    fir_filter_direct uut(
        .x(x),
{conns}
    );
    integer errors;
    initial begin
        errors = 0; x = 1; #10;
{checks}
        if (errors == 0) $display("ALL TESTS PASSED");
        else $display("FAILED: %0d errors", errors);
        $finish;
    end
endmodule
"""


def check_coefficients(payload: dict) -> tuple[list[int], int]:
    csd_list = [c["csd"] for c in payload["coefficients"]]
    raws, scale_power = align(csd_list)
    expected = []
    for coeff, raw in zip(payload["coefficients"], raws):
        assert coeff["value"] == parse_csd(coeff["csd"]), (
            f"tap {coeff['index']}: value != csd"
        )
        want = round(coeff["value"] * 2**scale_power)
        assert raw_to_int(raw) == want, f"tap {coeff['index']}: fixed-point misalignment"
        expected.append(want)
    return expected, scale_power


def run_case(cli: Path, spec: Path, work: Path, form: str) -> None:
    out_json = work / f"{form}.json"
    dut_v = work / f"{form}.v"
    subprocess.run(
        [str(cli), str(spec), "-o", str(out_json), "--verilog-out", str(dut_v)],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(out_json.read_text(encoding="utf-8"))
    expected, _ = check_coefficients(payload)
    verilog = payload["verilog"]

    if "posedge clk" in verilog:
        width = output_width(verilog, " y")
        tb = gen_transpose_tb(payload, expected, width)
    else:
        widths = [output_width(verilog, f" h{i}") for i in range(len(expected))]
        tb = gen_direct_tb(payload, expected, widths)

    tb_file = work / f"tb_{form}.v"
    tb_file.write_text(tb, encoding="utf-8")
    vvp_file = work / f"{form}.vvp"
    subprocess.run(
        [IVERILOG, "-g2012", "-o", str(vvp_file), str(tb_file), str(dut_v)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = subprocess.run(
        [VVP, str(vvp_file)], check=True, capture_output=True, text=True
    )
    combined = result.stdout + result.stderr
    if "ALL TESTS PASSED" not in combined or "FAIL" in combined:
        raise AssertionError(f"{form} simulation failed:\n{combined}")
    print(f"[PASS] {form}: {len(expected)} taps simulated")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", default=None, help="path to FirDesign executable")
    args = parser.parse_args()

    if IVERILOG is None or VVP is None:
        print("error: iverilog/vvp not installed", file=sys.stderr)
        sys.exit(1)

    cli = Path(args.cli) if args.cli else next(
        (
            p
            for p in (
                ROOT / "build" / "windows" / "x64" / "release" / "FirDesign.exe",
                ROOT / "build" / "linux" / "x86_64" / "release" / "FirDesign",
                ROOT / "build" / "standalone" / "Release" / "FirDesign.exe",
                ROOT / "build" / "standalone" / "FirDesign",
                ROOT / "build" / "Release" / "FirDesign.exe",
                ROOT / "build" / "FirDesign.exe",
                ROOT / "build" / "FirDesign",
            )
            if p.is_file()
        ),
        None,
    )
    if cli is None:
        print("error: FirDesign not found; pass --cli", file=sys.stderr)
        sys.exit(1)

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        run_case(cli, ROOT / "sample_filter.json", work, "transpose")
        run_case(cli, ROOT / "sample_filter_direct.json", work, "direct")


if __name__ == "__main__":
    main()
