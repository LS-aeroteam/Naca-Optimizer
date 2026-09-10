from scipy.optimize import minimize, differential_evolution
import numpy as np
import logging

from .airfoil import naca4_airfoil
from .panel_method import run_panel_analysis

# Configure basic logging
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

class NacaOptimizer:
    """
    Optimizes a NACA 4-digit airfoil to meet specific aerodynamic targets.
    """

    def __init__(self, reynolds, alpha, target_cl, max_height_box=0.3, chord=1.0, num_panels=160):
        self.reynolds = reynolds
        self.alpha = alpha
        self.target_cl = target_cl
        
        self.max_height_box = max_height_box
        self.chord = chord
        self.num_panels = num_panels
        
        self.eval_count = 0
        self.history = []

        # Penalty constants
        self.CL_ERROR_WEIGHT = 10.0
        self.BB_PENALTY_WEIGHT = 1e4
        self.FAILURE_PENALTY = 1e6

    def _objective_function(self, params):
        """
        The objective function to be minimized.
        Calculates a score based on potential Cl and bounding box constraints.
        """
        self.eval_count += 1
        m, p, t = params

        # Unpack parameters and check bounds
        if not (0.0 <= m <= 0.20 and 0.1 <= p <= 0.8 and 0.05 <= t <= 0.35):
            self.history.append([self.eval_count, m, p, t, "Bounds", "Bounds", self.FAILURE_PENALTY])
            return self.FAILURE_PENALTY

        # Generate airfoil with consistent number of points based on panels
        num_points = int(self.num_panels / 2) + 1
        X, Y, _ = naca4_airfoil(m, p, t, num_points=num_points)
        
        # Calculate bounding box physical dimensions (only height/thickness is relevant)
        height = (np.max(Y) - np.min(Y)) * self.chord
        
        bb_penalty = 0.0
        bb_flag = ""
        if height > self.max_height_box:
            bb_penalty += ((height - self.max_height_box) * self.BB_PENALTY_WEIGHT) ** 2
            bb_flag = "H"
        
        try:
            # Run panel analysis with verbose=False to avoid spamming the console
            panel_res = run_panel_analysis(X, Y, self.alpha, verbose=False)
            cl = panel_res['cl_potential']
        except Exception:
            cl = None

        if cl is not None:
            cl_error = ((cl - self.target_cl) * self.CL_ERROR_WEIGHT) ** 2
            score = cl_error + bb_penalty
            
            print(f"| {self.eval_count:4d} | {m:.4f} | {p:.4f} | {t:.4f} | {cl:.4f} | {bb_flag:^6} | {score:.4e} |")
            self.history.append([self.eval_count, m, p, t, cl, 0.0, score]) # Store 0.0 for Cd
            return score
        else:
            thickness_penalty = ((0.12 - t)**2 * 1e5)
            camber_penalty = ((0.05 - m)**2 * 1e5)
            emergency_score = self.FAILURE_PENALTY + bb_penalty + thickness_penalty + camber_penalty
            
            print(f"| {self.eval_count:4d} | {m:.4f} | {p:.4f} | {t:.4f} | {'Failed':^6} | {bb_flag:^6} | {emergency_score:.4e} |")
            self.history.append([self.eval_count, m, p, t, "Failed", "Failed", emergency_score])
            return emergency_score

    def optimize(self, initial_guess, bounds, max_iter=50):
        """
        Runs the hybrid optimization process (GA + SLSQP).

        Args:
            initial_guess (list): Initial values for [m, p, t].
            bounds (list of tuples): Bounds for [(m), (p), (t)].
            max_iter (int): Maximum number of iterations for the optimizer.

        Returns:
            scipy.optimize.OptimizeResult: The result object from scipy.minimize.
        """
        print("[+] Starting Hybrid Panel-Method Airfoil Optimization...")
        print("-" * 70)
        print(f"| {'Eval':^4} | {'m':^6} | {'p':^6} | {'t':^6} | {'Cl':^6} | {'BB OOB':^6} | {'Score':^10} |")
        print("-" * 70)

        self.eval_count = 0
        self.history = []

        # Phase 1: Global Search
        print("\n--> Phase 1: Genetic Algorithm (Global Search)")
        ga_result = differential_evolution(
            self._objective_function,
            bounds,
            maxiter=5, # Limit max iterations so it doesn't take forever
            popsize=5,
            polish=False # We will do our own polish with SLSQP
        )
        
        best_guess = ga_result.x
        
        # Phase 2: Local Refinement
        print("\n--> Phase 2: SLSQP (Local Refinement)")
        result = minimize(
            self._objective_function,
            best_guess,
            method='SLSQP',
            bounds=bounds,
            options={'disp': False, 'maxiter': max_iter, 'ftol': 1e-4, 'eps': 1e-4}
        )
        
        print("-" * 70)
        return result

    def get_optimization_history(self):
        return self.history
