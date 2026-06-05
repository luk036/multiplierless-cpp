"""Verify FIR filter coefficients meet passband/stopband specs."""
import json
import sys
import numpy as np


def verify_filter(json_path: str) -> bool:
    with open(json_path) as f:
        data = json.load(f)

    h = np.array([c["value"] for c in data["coefficients"]])
    N = len(h)

    # These should match whatever the user put in their JSON
    # Default to the standard values used in the example
    wpass_norm = 0.12
    wstop_norm = 0.20
    ripple = 0.125
    attn = 0.125

    delta1 = 20 * np.log10(1 + ripple)
    delta2 = 20 * np.log10(attn)
    L = pow(10, -delta1 / 20)
    U = pow(10, +delta1 / 20)
    S = pow(10, +delta2 / 20)
    Lsq = L * L
    Usq = U * U
    Ssq = S * S

    w = np.linspace(0, np.pi, 8192)
    H = np.zeros(len(w), dtype=complex)
    for i, wi in enumerate(w):
        H[i] = sum(h[k] * np.exp(-1j * k * wi) for k in range(N))
    R = np.abs(H) ** 2

    wpass_idx = int(wpass_norm * len(w))
    wstop_idx = int(wstop_norm * len(w))

    passband_min = np.min(R[:wpass_idx])
    passband_max = np.max(R[:wpass_idx])
    stopband_max = np.max(R[wstop_idx:])

    all_ok = True

    print(f"Filter: {N} taps, {data['csd_nnz']} CSD digits, {data['iterations']} iterations\n")
    print(f"Spec: passband 0-{wpass_norm}*pi, stopband {wstop_norm}*pi-pi")
    print(f"Bounds: L^2={Lsq:.4f}, U^2={Usq:.4f}, S^2={Ssq:.6f}\n")

    # Passband check
    if passband_min >= Lsq and passband_max <= Usq:
        print(f"[PASS] Passband: R(ω) ∈ [{passband_min:.4f}, {passband_max:.4f}]")
    else:
        print(f"[FAIL] Passband: R(ω) ∈ [{passband_min:.4f}, {passband_max:.4f}]")
        print(f"       Expected: [{Lsq:.4f}, {Usq:.4f}]")
        all_ok = False

    # Stopband check
    if stopband_max <= Ssq:
        print(f"[PASS] Stopband: max R(ω) = {stopband_max:.6f} ≤ {Ssq:.6f}")
    else:
        print(f"[FAIL] Stopband: max R(ω) = {stopband_max:.6f} > {Ssq:.6f}")
        all_ok = False

    # CSD digit check
    bad_taps = []
    nnz_target = data["csd_nnz"]
    for c in data["coefficients"]:
        nnz_count = c["csd"].count("+") + c["csd"].count("-")
        if nnz_count > nnz_target:
            bad_taps.append((c["index"], nnz_count))

    if not bad_taps:
        print(f"[PASS] CSD digits: all ≤ {nnz_target} non-zero per tap")
    else:
        print(f"[FAIL] CSD digits: {len(bad_taps)} taps exceed nnz={nnz_target}")
        for idx, count in bad_taps[:5]:
            print(f"       h[{idx}] has {count} non-zero digits")
        all_ok = False

    # Coefficient validity
    if np.all(np.isfinite(h)):
        print("[PASS] All coefficients finite")
    else:
        print("[FAIL] NaN/Inf in coefficients")
        all_ok = False

    print(f"\n{'ALL CHECKS PASSED' if all_ok else 'SOME CHECKS FAILED'}")
    return all_ok


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python verify_filter.py <output.json>")
        sys.exit(1)
    ok = verify_filter(sys.argv[1])
    sys.exit(0 if ok else 1)
