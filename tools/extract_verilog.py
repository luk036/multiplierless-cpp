"""Extract Verilog from FirDesign JSON output and save as standalone .v file.

Usage:
    python tools/extract_verilog.py                              # run default spec
    python tools/extract_verilog.py my_filter.json               # run custom spec
    python tools/extract_verilog.py output.json                  # extract from existing output
    python tools/extract_verilog.py output.json -o fir.v         # extract to named file
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path


def find_cli() -> Path:
    """Locate the FirDesign executable under the cmake build tree."""
    root = Path(__file__).resolve().parent.parent
    # Search common build output locations
    candidates = [
        root / "build" / "standalone" / "Debug" / "FirDesign.exe",
        root / "build" / "standalone" / "Release" / "FirDesign.exe",
        root / "build" / "windows" / "x64" / "release" / "FirDesign.exe",
    ]
    for p in candidates:
        if p.is_file():
            return p
    print("error: FirDesign.exe not found. Build the project first.", file=sys.stderr)
    sys.exit(1)


def extract(data: dict, output: Path) -> None:
    """Write the verilog field from JSON data to *output*."""
    verilog = data.get("verilog")
    if not verilog:
        print("error: no 'verilog' field in JSON data", file=sys.stderr)
        sys.exit(1)
    output.write_text(verilog, encoding="utf-8")
    iters = data.get("iterations", "?")
    print(f"Verilog → {output} ({len(verilog)} chars, "
          f"{verilog.count(chr(10))} lines, {iters} iterations)")


def run_cli(cli: Path, spec: Path) -> dict:
    """Run FirDesign on *spec* and return parsed JSON output."""
    result = subprocess.run([str(cli), str(spec)], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"error: FirDesign failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError as e:
        print(f"error: invalid JSON from FirDesign:\n{e}", file=sys.stderr)
        sys.exit(1)


def is_spec_file(path: Path) -> bool:
    """Heuristic: a spec file has 'filter_order' but NOT 'coefficients'."""
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        return isinstance(data, dict) and "filter_order" in data and "coefficients" not in data
    except (json.JSONDecodeError, OSError):
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", default="sample_filter.json",
                        help="Filter spec JSON or existing output JSON")
    parser.add_argument("-o", "--output", default=None,
                        help="Output .v file (default: <module_name>.v or fir_filter.v)")
    parser.add_argument("--cli", default=None,
                        help="Path to FirDesign executable (auto-detected if omitted)")
    args = parser.parse_args()

    input_path = Path(args.input)

    # Determine whether input is a spec (run CLI) or existing output (extract directly)
    if input_path.is_file() and is_spec_file(input_path):
        cli = Path(args.cli) if args.cli else find_cli()
        print(f"Running: {cli} {input_path}")
        data = run_cli(cli, input_path)
    elif input_path.is_file():
        data = json.loads(input_path.read_text(encoding="utf-8"))
    else:
        print(f"error: {input_path} not found", file=sys.stderr)
        sys.exit(1)

    # Determine output filename
    if args.output:
        output_path = Path(args.output)
    else:
        module = data.get("verilog", "")
        # Extract module name from the first line of Verilog
        name = "fir_filter"
        for line in module.splitlines():
            if line.startswith("module "):
                parts = line.split()
                if len(parts) >= 2:
                    name = parts[1]
                break
        output_path = Path(f"{name}.v")

    extract(data, output_path)

    # Also save the full JSON side-by-side (distinct name to avoid overwriting input)
    json_path = output_path.with_stem(output_path.stem + "_output").with_suffix(".json")
    json_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    print(f"JSON  → {json_path}")


if __name__ == "__main__":
    main()
