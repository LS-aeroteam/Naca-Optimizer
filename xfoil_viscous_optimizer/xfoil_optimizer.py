from naca_core.optimization import BaseNacaOptimizer
from naca_core.xfoil import XFoilAnalysis

# Largest accepted difference between the angle XFOIL converged at and the target angle (deg)
ALPHA_TOLERANCE = 0.1


class NacaOptimizer(BaseNacaOptimizer):
    """
    Finds a NACA 4-digit airfoil that reaches the target Cl with XFOIL (viscous).
    """
    SOLVER_NAME = "XFOIL"
    # XFOIL prints Cl with 4 decimals: with a step of 1e-4 the change of Cl due to p and t
    # is below 0.0001 and the finite-difference gradient becomes zero.
    FD_STEP = 2e-3

    CL_ERROR_WEIGHT = 10.0
    CD_PENALTY_WEIGHT = 1000.0

    def __init__(self, reynolds, alpha, target_cl, max_cd, mach=0.0, ncrit=9.0, chord=1.0,
                 max_height_box=None, num_panels=160, seed=None):
        super().__init__(chord=chord, max_height_box=max_height_box, seed=seed)
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
        cl_error = ((cl - self.target_cl) * self.CL_ERROR_WEIGHT) ** 2
        cd_penalty = (max(0, cd - self.max_cd) * self.CD_PENALTY_WEIGHT) ** 2
        return cl_error + cd_penalty
