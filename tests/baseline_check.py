"""
Numerical baseline for the in-house solver (geometry + panel method).

Run from anywhere:
    python tests/baseline_check.py --generate   # save the reference from the current code (once)
    python tests/baseline_check.py              # compare the current code with the saved reference

Exit code: 0 = every value matches, 1 = at least one mismatch or missing reference.
"""
import argparse
import json
import os
import sys

import numpy as np

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)
sys.path.insert(0, os.path.join(REPO_ROOT, "inhouse_potential_optimizer"))

# --- Modules under test: when files move (refactor), update only these two lines ---
from inhouse_core.airfoil import naca4_airfoil
from inhouse_core.panel_method import run_panel_analysis
# ------------------------------------------------------------------------------------

REFERENCE_FILE = os.path.join(TESTS_DIR, "baseline_reference.json")
TOLERANCE = 1e-9          # max absolute difference allowed for every value
NUM_PANELS = 160
AIRFOILS = {"0012": (0.00, 0.0, 0.12), "2412": (0.02, 0.4, 0.12), "4412": (0.04, 0.4, 0.12)}
ALPHAS = [-4.0, 0.0, 5.0]


def compute_current():
    """Runs geometry and panel method on the test matrix and returns plain lists/floats."""
    data = {}
    for name, (m, p, t) in AIRFOILS.items():
        X, Y, _ = naca4_airfoil(m, p, t, num_points=NUM_PANELS // 2 + 1)
        data[f"NACA {name} | X"] = X.tolist()
        data[f"NACA {name} | Y"] = Y.tolist()
        for alpha in ALPHAS:
            res = run_panel_analysis(X, Y, alpha, verbose=False)
            data[f"NACA {name} | alpha {alpha:+.1f} | Cl"] = float(res["cl_potential"])
            data[f"NACA {name} | alpha {alpha:+.1f} | Cp"] = np.asarray(res["Cp"]).tolist()
    return data


def generate():
    if os.path.exists(REFERENCE_FILE):
        print(f"[!] Reference already exists: {REFERENCE_FILE}")
        print("    It is not overwritten, so a refactor cannot silently change the baseline.")
        print("    Delete the file by hand only if the change in results is intended.")
        return 1
    with open(REFERENCE_FILE, "w") as f:
        json.dump(compute_current(), f, indent=1)
    print(f"[+] Reference saved to: {REFERENCE_FILE}")
    return 0


def check():
    if not os.path.exists(REFERENCE_FILE):
        print("[!] No reference found. Run first: python tests/baseline_check.py --generate")
        return 1
    with open(REFERENCE_FILE) as f:
        reference = json.load(f)
    current = compute_current()

    print("=" * 70)
    print(f"{'Quantity':<40} | {'Max |diff|':>11} | Result")
    print("-" * 70)
    failures = 0
    for key in reference:
        if key not in current:
            print(f"{key:<40} | {'missing':>11} | FAIL")
            failures += 1
            continue
        ref = np.asarray(reference[key], dtype=float)
        cur = np.asarray(current[key], dtype=float)
        if ref.shape != cur.shape:
            print(f"{key:<40} | {'shape':>11} | FAIL")
            failures += 1
            continue
        diff = float(np.max(np.abs(ref - cur)))
        ok = diff <= TOLERANCE
        failures += not ok
        print(f"{key:<40} | {diff:11.2e} | {'ok' if ok else 'FAIL'}")
    print("-" * 70)
    if failures:
        print(f"[!] {failures} quantities differ from the baseline (tolerance {TOLERANCE:g}).")
        return 1
    print(f"[+] All {len(reference)} quantities match the baseline (tolerance {TOLERANCE:g}).")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Numerical baseline check for the in-house solver.")
    parser.add_argument("--generate", action="store_true", help="save the reference from the current code")
    args = parser.parse_args()
    sys.exit(generate() if args.generate else check())
