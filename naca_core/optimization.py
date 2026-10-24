"""
Shared optimization logic for both solvers.

Both optimizers use the same two-phase search (Genetic Algorithm = scipy
differential evolution, then SLSQP), the same bounding-box rule and the same
terminal table. Each solver only defines how an airfoil is analysed and how the
score is computed. The airfoil family (naca_core.airfoil.AirfoilFamily) gives the
parameters, their search limits, the starting airfoil and the geometry.
"""
import numpy as np
from scipy.optimize import minimize, differential_evolution

from .airfoil import NACA4

# Score given to failed analyses and to airfoils outside the bounding box. It must be larger
# than the score of any valid airfoil, including the objective penalties of the XFOIL solver
# (up to about 4e8 for a Cl 2.0 away from the target), so a valid airfoil always ranks better.
FAILURE_PENALTY = 1e9
BOX_PENALTY_WEIGHT = 1000.0

# Terminal table: one definition for both solvers, one column per parameter of the family
def _table_header(params, width=6):
    columns = "".join(f"{name:^{width}} | " for name in params)
    return f"| {'Eval':^4} | {columns}{'Cl':^7} | {'Cd':^7} | {'BB':^4} | {'Score':^11} |"


def print_phase(title):
    """Prints a phase header centred in a 70-character frame."""
    print("\n" + "=" * 70)
    print(title.center(70).rstrip())
    print("=" * 70)


class BaseNacaOptimizer:
    """
    Common part of the in-house and XFOIL optimizers.

    Subclasses implement:
        _analyse(X, Y) -> (cl, cd, status)   cl/cd are None on failure,
                                             status is a short label ("Failed", "Timeout")
        _score(cl, cd) -> float              score to minimise for a successful analysis
    and set SOLVER_NAME, FD_STEP (finite-difference step used by SLSQP) and NUM_POINTS.

    History rows: [eval, *parameters of the family, Cl, Cd, score].
    """
    SOLVER_NAME = "Solver"
    FD_STEP = 1e-4
    NUM_POINTS = 100  # points per surface used to build the airfoil for the analysis

    def __init__(self, chord=1.0, max_height_box=None, seed=None, family=NACA4):
        self.chord = chord
        self.max_height_box = max_height_box
        self.seed = seed
        self.family = family
        # Parameter columns: 6 characters (0.0200), 7 when a parameter can be negative (-0.1504)
        self.col_width = 7 if any(lo < 0 for lo, _ in family.bounds) else 6
        self.eval_count = 0
        self.history = []

    # --- to be implemented by each solver ---
    def _analyse(self, X, Y):
        raise NotImplementedError

    def _score(self, cl, cd):
        raise NotImplementedError

    # --- shared logic ---
    def _print_row(self, n, values, cl, cd, bb, score):
        cl_str = cl if isinstance(cl, str) else f"{cl:.4f}"
        cd_str = cd if isinstance(cd, str) else f"{cd:.5f}"
        columns = "".join(f"{v:{self.col_width}.4f} | " for v in values)
        print(f"| {n:4d} | {columns}{cl_str:>7} | {cd_str:>7} | {bb:^4} | {score:11.4e} |")

    def _objective_function(self, params):
        self.eval_count += 1
        # Keep parameters inside the limits (guards against round-off at the bounds)
        values = [float(np.clip(v, lo, hi)) for v, (lo, hi) in zip(params, self.family.bounds)]

        # Parameters that give no real airfoil (e.g. CST upper surface below the lower one): not analysed
        shape_error = self.family.shape_error(values)
        if shape_error:
            score = self.family.failure_score(values, FAILURE_PENALTY)
            self._print_row(self.eval_count, values, shape_error, "-", "", score)
            self.history.append([self.eval_count, *values, shape_error, shape_error, score])
            return score

        X, Y = self.family.coordinates(values, self.NUM_POINTS)

        # Bounding box: geometric check only, the solver is not run if the airfoil is too tall
        if self.max_height_box is not None:
            height = (np.max(Y) - np.min(Y)) * self.chord
            if height > self.max_height_box:
                score = FAILURE_PENALTY + ((height - self.max_height_box) * BOX_PENALTY_WEIGHT) ** 2
                self._print_row(self.eval_count, values, "-", "-", "OUT", score)
                self.history.append([self.eval_count, *values, "Box Fail", "Box Fail", score])
                return score

        cl, cd, status = self._analyse(X, Y)

        if cl is None:
            # Failed analysis: large penalty, slightly smaller for more usual shapes
            score = self.family.failure_score(values, FAILURE_PENALTY)
            self._print_row(self.eval_count, values, status, "-", "", score)
            self.history.append([self.eval_count, *values, status, status, score])
            return score

        score = self._score(cl, cd)
        self._print_row(self.eval_count, values, cl, "-" if cd is None else cd, "", score)
        self.history.append([self.eval_count, *values, cl, "" if cd is None else cd, score])
        return score

    def optimize(self, max_iter=50):
        """
        Two-phase search: Genetic Algorithm (differential evolution, global) then SLSQP (local).

        Returns:
            scipy.optimize.OptimizeResult: result of the SLSQP phase.
        """
        header = _table_header(self.family.params, self.col_width)
        separator = "-" * len(header)
        print(f"[+] Starting Hybrid {self.SOLVER_NAME} Airfoil Optimization...")
        print(separator)
        print(header)
        print(separator)

        self.eval_count = 0
        self.history = []

        print("\n--> Phase 1: Genetic Algorithm (Global Search)")
        ga_result = differential_evolution(
            self._objective_function,
            self.family.bounds,
            maxiter=5,       # short global search: the local phase refines the result
            popsize=5,
            polish=False,    # the local refinement is done below with SLSQP
            x0=self.family.initial_guess,
            rng=np.random.default_rng(self.seed),
        )

        print("\n--> Phase 2: SLSQP (Local Refinement)")
        result = minimize(
            self._objective_function,
            ga_result.x,
            method='SLSQP',
            bounds=self.family.bounds,
            options={'disp': False, 'maxiter': max_iter, 'ftol': 1e-4, 'eps': self.FD_STEP}
        )

        print(separator)
        return result

    def get_optimization_history(self):
        return self.history
