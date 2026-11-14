"""
Fixed case of the terminal non-regression test (tests/terminal_regression_check.py).

Every input is written out, so the test does not depend on the defaults in naca_core/input_file.py.
Do not change these values: the saved reference output was produced with them.
In-house solver only: XFOIL output depends on the XFOIL executable and the platform.
"""

SOLVER = "inhouse"

FLUID = "air"
SPEED = 50.0
CHORD = 1.0
ALPHA = 4.0

OBJECTIVE = "min_cd"
TARGET_CL = 0.8
MAX_CD = 0.02

MAX_HEIGHT_BOX = 0.3
NUM_PANELS = 160
SEED = 2412

SPAN = 1.0
