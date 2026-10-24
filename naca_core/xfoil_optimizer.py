from .airfoil import NACA4
from .optimization import BaseNacaOptimizer
from .xfoil import XFoilAnalysis

# Largest accepted difference between the angle XFOIL converged at and the target angle (deg)
ALPHA_TOLERANCE = 0.1

# Objectives
MIN_CD = "min_cd"   # minimum Cd with Cl = target Cl +/- CL_TOLERANCE
MAX_CL = "max_cl"   # maximum Cl with Cd <= max Cd
CL_TOLERANCE = 0.005


class NacaOptimizer(BaseNacaOptimizer):
    """
    Optimizes an airfoil of the given family (default NACA 4-digit) with XFOIL (viscous) for one of two objectives:
        MIN_CD: minimum Cd with Cl inside target_cl +/- CL_TOLERANCE
        MAX_CL: maximum Cl with Cd not above max_cd
    """
    SOLVER_NAME = "XFOIL"
    # XFOIL prints Cl with 4 decimals: with a step of 1e-4 the change of Cl due to p and t
    # is below 0.0001 and the finite-difference gradient becomes zero.
    FD_STEP = 2e-3

    # Penalty weights. Scores are in drag counts (1 count = 0.0001 of Cd) for MIN_CD.
    # Linear term: its slope is already there at the edge of the tolerance, where the quadratic
    # term is flat, so leaving the tolerance never pays off: 0.0001 of Cl outside it costs
    # 100 counts, against about 3 counts of Cd saved (seed 657086: 0.0003 of Cl outside the
    # tolerance saved 9 counts of Cd, while the quadratic term charged only 9 counts for it).
    CL_LINEAR_WEIGHT = 1e6
    # Quadratic term: 0.001 of Cl outside the tolerance costs 100 counts more, growing fast further out.
    CL_PENALTY_WEIGHT = 1e4
    # MAX_CL, same idea. Linear term: 1 count of Cd above the limit costs 1.0 of Cl from the edge
    # of the limit, against about 0.024 of Cl gained per count near it (default conditions, Cd 0.006).
    CD_LINEAR_WEIGHT = 1e4
    # Quadratic term: 1 count of Cd above the limit costs 1.0 of Cl more, growing fast further out.
    CD_PENALTY_WEIGHT = 1e4

    def __init__(self, reynolds, alpha, objective=MIN_CD, target_cl=None, max_cd=None, mach=0.0, ncrit=9.0,
                 chord=1.0, max_height_box=None, num_panels=160, seed=None, family=NACA4):
        super().__init__(chord=chord, max_height_box=max_height_box, seed=seed, family=family)
        if objective == MIN_CD and target_cl is None:
            raise ValueError("MIN_CD objective needs target_cl")
        if objective == MAX_CL and max_cd is None:
            raise ValueError("MAX_CL objective needs max_cd")
        self.objective = objective
        self.reynolds = reynolds
        self.alpha = alpha
        self.target_cl = target_cl
        self.max_cd = max_cd
        self.mach = mach
        self.ncrit = ncrit
        self.num_panels = num_panels
        self.NUM_POINTS = int(num_panels / 2) + 1  # XFOIL analyses these points directly (no PANE)

    def _analyse(self, X, Y):
        analysis = XFoilAnalysis(
            airfoil_name=f"opt_naca_{self.eval_count}",
            alpha=self.alpha,
            reynolds=self.reynolds,
            mach=self.mach,
            ncrit=self.ncrit,
        )
        cl, cd, achieved_alpha = analysis.run_analysis(X, Y)

        if analysis.last_error == "timeout":
            return None, None, "Timeout"
        if cl is None or cd is None or cd <= 0:
            return None, None, "Failed"
        # A Cl computed at another angle is not a valid result: treated as a failed analysis
        if abs(achieved_alpha - self.alpha) > ALPHA_TOLERANCE:
            return None, None, "Failed"
        return cl, cd, None

    def _score(self, cl, cd):
        if self.objective == MIN_CD:
            # Inside the Cl tolerance only Cd counts; outside, a linear + quadratic penalty is added
            cl_excess = max(0.0, abs(cl - self.target_cl) - CL_TOLERANCE)
            return cd * 1e4 + cl_excess * self.CL_LINEAR_WEIGHT + (cl_excess * self.CL_PENALTY_WEIGHT) ** 2
        # MAX_CL: the higher Cl, the lower the score; Cd above the limit gets a linear + quadratic penalty
        cd_excess = max(0.0, cd - self.max_cd)
        return -cl + cd_excess * self.CD_LINEAR_WEIGHT + (cd_excess * self.CD_PENALTY_WEIGHT) ** 2
