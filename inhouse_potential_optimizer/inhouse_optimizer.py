from naca_core.optimization import BaseNacaOptimizer
from naca_core.panel_method import run_panel_analysis


class NacaOptimizer(BaseNacaOptimizer):
    """
    Finds a NACA 4-digit airfoil that reaches the target Cl with the in-house panel method.
    Potential flow: Cd is not available.
    """
    SOLVER_NAME = "Panel-Method"
    FD_STEP = 1e-4  # the panel method is smooth: a small finite-difference step is accurate

    CL_ERROR_WEIGHT = 10.0

    def __init__(self, reynolds, alpha, target_cl, max_height_box=0.3, chord=1.0, num_panels=160, seed=None):
        super().__init__(chord=chord, max_height_box=max_height_box, seed=seed)
        self.reynolds = reynolds
        self.alpha = alpha
        self.target_cl = target_cl
        self.num_panels = num_panels
        self.NUM_POINTS = int(num_panels / 2) + 1  # same paneling as the final analysis

    def _analyse(self, X, Y):
        try:
            res = run_panel_analysis(X, Y, self.alpha, verbose=False)
            return res['cl_potential'], None, None
        except Exception:
            return None, None, "Failed"

    def _score(self, cl, cd):
        return ((cl - self.target_cl) * self.CL_ERROR_WEIGHT) ** 2
