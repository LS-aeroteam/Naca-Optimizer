"""
NACA Optimizer - input file
===========================

    python run.py --input inputs.py

How to use it
-------------
1. Change the values below and save.
2. Run the command above: no question is asked in the terminal.

- Every line is optional: delete it (or put a # in front of it) and the default is used.
- Text goes "in quotes", numbers without. None means "not set".
- Keep several cases by copying this file (e.g. cases/water_6deg.py) and passing the copy to --input.
- Defaults are the values written below; units are in [brackets].
"""

# =============================================================================
#  1. SOLVER
# =============================================================================
#   "inhouse"  our panel method: potential flow, fast, no drag (Cd)
#   "xfoil"    XFOIL: viscous, slower, gives drag (Cd)
#   The command-line option --solver, if given, replaces this value.

SOLVER = "inhouse"


# =============================================================================
#  2. FLOW CONDITIONS          Reynolds and Mach are computed from these values
# =============================================================================

FLUID = "air"           # "air"    sea level, nu = 1.46e-5 m^2/s, a = 340.3 m/s
                        # "water"  20 deg C,  nu = 1.00e-6 m^2/s, a = 1482 m/s
SPEED = 50.0            # [m/s]  design speed
CHORD = 1.0             # [m]    airfoil chord
ALPHA = 4.0             # [deg]  angle of attack


# =============================================================================
#  3. OBJECTIVE
# =============================================================================
#   OBJECTIVE is used only with SOLVER = "xfoil".
#   With "inhouse" it is ignored: the in-house solver always reaches TARGET_CL.

OBJECTIVE = "min_cd"    # "min_cd"  minimum Cd with Cl = TARGET_CL +/- 0.005
                        # "max_cl"  maximum Cl with Cd <= MAX_CD
TARGET_CL = 0.8         # [-]    target lift coefficient   (used by "min_cd" and "inhouse")
MAX_CD = 0.006          # [-]    maximum drag coefficient  (used by "max_cl")


# =============================================================================
#  4. SEARCH
# =============================================================================

MAX_HEIGHT_BOX = 0.3    # [m]    maximum total height of the airfoil (thickness + camber)
NUM_PANELS = 160        # [-]    60 = fast, 100 = medium, 160 = accurate
SEED = None             # [-]    None = random (printed at the end of the run);
                        #        a whole number (e.g. 42) repeats a run exactly


# =============================================================================
#  5. EXPORT
# =============================================================================

SPAN = None             # [m]    span of the extruded wing in the STL file; None = same as CHORD
