"""
Shared optimization logic for both solvers.

Both optimizers use the same search limits, the same two-phase search
(Genetic Algorithm = scipy differential evolution, then SLSQP), the same
bounding-box rule and the same terminal table. Each solver only defines how
an airfoil is analysed and how the score is computed.
"""
import numpy as np
from scipy.optimize import minimize, differential_evolution

from .airfoil import naca4_airfoil

# Search limits for [m, p, t] and starting airfoil (NACA 2412), shared by both solvers
BOUNDS = [(0.0, 0.09), (0.1, 0.7), (0.05, 0.25)]
INITIAL_GUESS = [0.02, 0.4, 0.12]

# Score given to failed analyses and to airfoils outside the bounding box. It must be larger
# than the score of any valid airfoil, including the objective penalties of the XFOIL solver
# (up to about 4e8 for a Cl 2.0 away from the target), so a valid airfoil always ranks better.
FAILURE_PENALTY = 1e9
BOX_PENALTY_WEIGHT = 1000.0

# Terminal table: one definition for both solvers
_TABLE_HEADER = (f"| {'Eval':^4} | {'m':^6} | {'p':^6} | {'t':^6} | {'Cl':^7} | "
                 f"{'Cd':^7} | {'BB':^4} | {'Score':^11} |")
_SEPARATOR = "-" * len(_TABLE_HEADER)


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
    """
    SOLVER_NAME = "Solver"
    FD_STEP = 1e-4
    NUM_POINTS = 100  # points per surface used to build the airfoil for the analysis

    def __init__(self, chord=1.0, max_height_box=None, seed=None):
        self.chord = chord
        self.max_height_box = max_height_box
        self.seed = seed
        self.eval_count = 0
        self.history = []

    # --- to be implemented by each solver ---
    def _analyse(self, X, Y):
        raise NotImplementedError

    def _score(self, cl, cd):
        raise NotImplementedError

    # --- shared logic ---
    @staticmethod
    def _print_row(n, m, p, t, cl, cd, bb, score):
        cl_str = cl if isinstance(cl, str) else f"{cl:.4f}"
        cd_str = cd if isinstance(cd, str) else f"{cd:.5f}"
        print(f"| {n:4d} | {m:6.4f} | {p:6.4f} | {t:6.4f} | {cl_str:>7} | {cd_str:>7} | {bb:^4} | {score:11.4e} |")

    def _objective_function(self, params):
        self.eval_count += 1
        # Keep parameters inside the limits (guards against round-off at the bounds)
        m, p, t = (float(np.clip(v, lo, hi)) for v, (lo, hi) in zip(params, BOUNDS))

        X, Y, _ = naca4_airfoil(m, p, t, num_points=self.NUM_POINTS)

        # Bounding box: geometric check only, the solver is not run if the airfoil is too tall
        if self.max_height_box is not None:
            height = (np.max(Y) - np.min(Y)) * self.chord
            if height > self.max_height_box:
                score = FAILURE_PENALTY + ((height - self.max_height_box) * BOX_PENALTY_WEIGHT) ** 2
                self._print_row(self.eval_count, m, p, t, "-", "-", "OUT", score)
                self.history.append([self.eval_count, m, p, t, "Box Fail", "Box Fail", score])
                return score

        cl, cd, status = self._analyse(X, Y)

        if cl is None:
            # Failed analysis: large penalty, slightly smaller for more usual shapes
            score = FAILURE_PENALTY + (0.12 - t) ** 2 * 1e5 + (0.05 - m) ** 2 * 1e5
            self._print_row(self.eval_count, m, p, t, status, "-", "", score)
            self.history.append([self.eval_count, m, p, t, status, status, score])
            return score

        score = self._score(cl, cd)
        self._print_row(self.eval_count, m, p, t, cl, "-" if cd is None else cd, "", score)
        self.history.append([self.eval_count, m, p, t, cl, "" if cd is None else cd, score])
        return score

    def optimize(self, max_iter=50):
        """
        Two-phase search: Genetic Algorithm (differential evolution, global) then SLSQP (local).

        Returns:
            scipy.optimize.OptimizeResult: result of the SLSQP phase.
        """
        print(f"[+] Starting Hybrid {self.SOLVER_NAME} Airfoil Optimization...")
        print(_SEPARATOR)
        print(_TABLE_HEADER)
        print(_SEPARATOR)

        self.eval_count = 0
        self.history = []

        print("\n--> Phase 1: Genetic Algorithm (Global Search)")
        ga_result = differential_evolution(
            self._objective_function,
            BOUNDS,
            maxiter=5,       # short global search: the local phase refines the result
            popsize=5,
            polish=False,    # the local refinement is done below with SLSQP
            x0=INITIAL_GUESS,
            rng=np.random.default_rng(self.seed),
        )

        print("\n--> Phase 2: SLSQP (Local Refinement)")
        result = minimize(
            self._objective_function,
            ga_result.x,
            method='SLSQP',
            bounds=BOUNDS,
            options={'disp': False, 'maxiter': max_iter, 'ftol': 1e-4, 'eps': self.FD_STEP}
        )

        print(_SEPARATOR)
        return result

    def get_optimization_history(self):
        return self.history
