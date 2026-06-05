"""Extract Verilog from FirDesign JSON output and save as standalone .v file."""
import json
import subprocess

# --- C++ CLI ---
result = subprocess.run(
    [r".\build\windows\x64\release\FirDesign.exe", "sample_filter.json"],
    capture_output=True, text=True,
)
data_cpp = json.loads(result.stdout)
verilog_cpp = data_cpp["verilog"]

with open("fir_filter_cpp.v", "w", encoding="utf-8") as f:
    f.write(verilog_cpp)
print(f"C++ Verilog → fir_filter_cpp.v ({len(verilog_cpp)} chars, "
      f"{verilog_cpp.count(chr(10))} lines, "
      f"{data_cpp['iterations']} iterations)")

# Keep full JSON output
with open("fir_design_output_cpp.json", "w", encoding="utf-8") as f:
    json.dump(data_cpp, f, indent=2)

# --- Python CLI ---
result = subprocess.run(
    ["conda", "run", "-n", "pyflow", "fir-design", "sample_filter.json"],
    capture_output=True, text=True,
    cwd=r"D:\github\py\multiplierless",
)
data_py = json.loads(result.stdout)
verilog_py = data_py["verilog"]

with open("fir_filter_py.v", "w", encoding="utf-8") as f:
    f.write(verilog_py)
print(f"Python Verilog → fir_filter_py.v ({len(verilog_py)} chars, "
      f"{verilog_py.count(chr(10))} lines, "
      f"{data_py['iterations']} iterations)")

# Keep full JSON output
with open("fir_design_output_py.json", "w", encoding="utf-8") as f:
    json.dump(data_py, f, indent=2)

print("Done.")
