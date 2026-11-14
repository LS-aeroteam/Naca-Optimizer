"""
Terminal non-regression test: the whole terminal output of a fixed run must not change.

Runs, as a user would:
    python run.py --input tests/regression_inputs.py      (in-house solver, seed 2412)
and compares its output, line by line, with the saved reference. It covers what
baseline_check.py does not: the optimization (every evaluated airfoil, its Cl and score),
the final airfoil, the messages and the export.

Run from anywhere:
    python tests/terminal_regression_check.py --generate   # save the reference from the current code (once)
    python tests/terminal_regression_check.py              # compare the current code with the saved reference

The run writes its results in Results/Re3424658_Alpha4.0_Cl0.8_inhouse_seed2412/, like any run.
Exit code: 0 = identical output, 1 = at least one line differs, missing reference or failed run.
"""
import argparse
import difflib
import os
import subprocess
import sys

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)

INPUT_FILE = "tests/regression_inputs.py"  # relative to REPO_ROOT: it is printed in the output
REFERENCE_FILE = os.path.join(TESTS_DIR, "terminal_regression_reference.txt")
TIMEOUT = 300       # [s] the run takes a few seconds
MAX_DIFF_LINES = 60  # lines of the diff printed on a mismatch


def run_and_capture():
    """Runs run.py with the fixed input file and returns its terminal output, normalised."""
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    proc = subprocess.run(
        [sys.executable, "run.py", "--input", INPUT_FILE],
        cwd=REPO_ROOT, stdin=subprocess.DEVNULL, capture_output=True,
        text=True, encoding="utf-8", timeout=TIMEOUT, env=env,
    )
    output = proc.stdout + proc.stderr
    if proc.returncode != 0:
        print(output)
        raise RuntimeError(f"run.py ended with exit code {proc.returncode}")
    # The repository path and the path separator depend on the machine: replace them
    output = output.replace(REPO_ROOT, "<repo>").replace("\\", "/")
    return [line.rstrip() for line in output.splitlines()]


def generate():
    if os.path.exists(REFERENCE_FILE):
        print(f"[!] Reference already exists: {REFERENCE_FILE}")
        print("    It is not overwritten, so a refactor cannot silently change the reference.")
        print("    Delete the file by hand only if the change in the terminal output is intended.")
        return 1
    lines = run_and_capture()
    with open(REFERENCE_FILE, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[+] Reference saved to: {REFERENCE_FILE} ({len(lines)} lines)")
    return 0


def check():
    if not os.path.exists(REFERENCE_FILE):
        print("[!] No reference found. Run first: python tests/terminal_regression_check.py --generate")
        return 1
    with open(REFERENCE_FILE, encoding="utf-8") as f:
        reference = f.read().splitlines()
    print(f"[+] Running: python run.py --input {INPUT_FILE}")
    current = run_and_capture()

    if current == reference:
        print(f"[+] Terminal output identical to the reference ({len(reference)} lines).")
        return 0

    diff = list(difflib.unified_diff(reference, current, "reference", "current", lineterm="", n=1))
    changed = sum(1 for line in diff if line.startswith(("-", "+")) and not line.startswith(("---", "+++")))
    print(f"[!] Terminal output differs from the reference ({changed} lines removed/added):")
    for line in diff[:MAX_DIFF_LINES]:
        print("    " + line)
    if len(diff) > MAX_DIFF_LINES:
        print(f"    ... {len(diff) - MAX_DIFF_LINES} more diff lines")
    return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Terminal non-regression test (fixed in-house run).")
    parser.add_argument("--generate", action="store_true", help="save the reference from the current code")
    args = parser.parse_args()
    sys.exit(generate() if args.generate else check())
