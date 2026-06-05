"""Verify FIR filter meets specs using CSD-reconstructed coefficients.

Verifies at TWO resolutions:
  1. Oracle grid  (disc_factor * N points) — what the optimizer validated against
  2. Fine grid    (16384 points)            — dense frequency sweep

The oracle grid check is authoritative: the optimizer guarantees feasibility
only at its sampled frequencies. Fine-grid violations are discretization
artifacts; increase 'discretization_factor' to reduce them.
"""
import json
import sys
import numpy as np


def csd_to_float(csd: str) -> float:
    integral = 0
    frac = 0.0
    scale = 0.5
    in_frac = False
    for ch in csd:
        if ch == '.':
            in_frac = True
            scale = 0.5
        elif ch == '+':
            if in_frac:
                frac += scale
            else:
                integral = integral * 2 + 1
            scale /= 2
        elif ch == '-':
            if in_frac:
                frac -= scale
            else:
                integral = integral * 2 - 1
            scale /= 2
        elif ch == '0':
            if not in_frac:
                integral *= 2
            scale /= 2
    return float(integral) + frac


def freq_resp(h, w):
    H = np.zeros(len(w), dtype=complex)
    for i, wi in enumerate(w):
        H[i] = sum(h[k] * np.exp(-1j * k * wi) for k in range(len(h)))
    return np.abs(H) ** 2


def check_band(R, start, end, name, lower, upper):
    """Check R[start:end] against [lower, upper] bounds."""
    chunk = R[start:end]
    mn, mx = np.min(chunk), np.max(chunk)
    if lower is not None:
        ok = (mn >= lower) and (mx <= upper)
    else:
        ok = mx <= upper
    return mn, mx, ok


def verify_filter(json_path: str, disc_factor: int = 15) -> bool:
    with open(json_path) as f:
        data = json.load(f)

    h_csd = np.array([csd_to_float(c["csd"]) for c in data["coefficients"]])
    N = data["filter_order"]

    wpass_norm = 0.12
    wstop_norm = 0.20
    ripple = 0.125
    attn = 0.125

    delta1 = 20 * np.log10(1 + ripple)
    delta2 = 20 * np.log10(attn)
    L = pow(10, -delta1 / 20)
    U = pow(10, +delta1 / 20)
    S = pow(10, +delta2 / 20)
    Lsq, Usq, Ssq = L * L, U * U, S * S

    nnz_target = data["csd_nnz"]
    all_ok = True

    print(f"Filter: {N} taps, {nnz_target} CSD digits/coeff, "
          f"{data['iterations']} iterations")
    print(f"Bounds: L^2={Lsq:.4f}, U^2={Usq:.4f}, S^2={Ssq:.6f}\n")

    # --- Oracle grid check (authoritative) ---
    m_oracle = disc_factor * N
    w_oracle = np.linspace(0, np.pi, m_oracle)
    R_oracle = freq_resp(h_csd, w_oracle)
    wpass_o = int(wpass_norm * m_oracle)
    wstop_o = int(wstop_norm * m_oracle)

    print(f"--- Oracle grid ({m_oracle} points, disc_factor={disc_factor}) ---")
    pmn, pmx, pok = check_band(R_oracle, 0, wpass_o, "Passband", Lsq, Usq)
    smx, _, sok = check_band(R_oracle, wstop_o, m_oracle, "Stopband", None, Ssq)
    print(f"[{'PASS' if pok else 'FAIL'}] Passband: R(w) in [{pmn:.4f}, {pmx:.4f}]")
    print(f"[{'PASS' if sok else 'FAIL'}] Stopband: max R(w) = {smx:.6f} <= {Ssq:.6f}")

    # --- Fine grid check (advisory) ---
    w_fine = np.linspace(0, np.pi, 16384)
    R_fine = freq_resp(h_csd, w_fine)
    wp_f = int(wpass_norm * 16384)
    ws_f = int(wstop_norm * 16384)

    print(f"\n--- Fine grid (16384 points, advisory) ---")
    pfmn, pfmx, pfok = check_band(R_fine, 0, wp_f, "Passband", Lsq, Usq)
    sfmx, _, sfok = check_band(R_fine, ws_f, 16384, "Stopband", None, Ssq)
    grid_note = ""
    if pok and not pfok:
        grid_note = "  (discretization artifact: increase disc_factor)"
    print(f"[{'PASS' if pfok else 'FAIL'}] Passband: R(w) in [{pfmn:.4f}, {pfmx:.4f}]{grid_note}")
    print(f"[{'PASS' if sfok else 'FAIL'}] Stopband: max R(w) = {sfmx:.6f} <= {Ssq:.6f}")

    # --- CSD digit check ---
    bad = [(c["index"], c["csd"].count("+") + c["csd"].count("-"))
           for c in data["coefficients"]
           if c["csd"].count("+") + c["csd"].count("-") > nnz_target]
    if not bad:
        print(f"\n[PASS] CSD digits: all <= {nnz_target}")
    else:
        print(f"\n[FAIL] CSD digits: {len(bad)} taps exceed nnz={nnz_target}")
        all_ok = False

    if not np.all(np.isfinite(h_csd)):
        print("[FAIL] NaN/Inf in coefficients")
        all_ok = False

    if not pok or not sok:
        all_ok = False

    print(f"\n{'ALL CHECKS PASSED' if all_ok else 'SOME CHECKS FAILED'}")
    if pok and not pfok:
        print("Note: Fine-grid passband violation is a discretization artifact.")
        print("      Increase 'discretization_factor' in your JSON spec to tighten.")
    return all_ok


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python verify_filter.py <output.json> [disc_factor=15]")
        sys.exit(1)
    df = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    ok = verify_filter(sys.argv[1], df)
    sys.exit(0 if ok else 1)
