from scipy.optimize import minimize
import numpy as np
import logging

from .airfoil import naca4_airfoil
from .xfoil import XFoilAnalysis

# Configure basic logging
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

class NacaOptimizer:
    """
    Optimizes a NACA 4-digit airfoil to meet specific aerodynamic targets.
    """

    def __init__(self, reynolds, alpha, target_cl, max_cd, mach=0.0, ncrit=9.0, chord=1.0, max_height=None):
        self.reynolds = reynolds
        self.alpha = alpha
        self.target_cl = target_cl
        self.max_cd = max_cd
        self.mach = mach
        self.ncrit = ncrit
        self.chord = chord
        self.max_height = max_height
        
        self.eval_count = 0
        self.history = []

        # Penalty constants
        self.CL_ERROR_WEIGHT = 10.0
        self.CD_PENALTY_WEIGHT = 1000.0
        self.ALPHA_PENALTY_WEIGHT = 1000.0
        self.FAILURE_PENALTY = 1e6

    def _objective_function(self, params):
        """
        The objective function to be minimized by SLSQP.
        Calculates a score based on aerodynamic performance.
        """
        self.eval_count += 1
        m, p, t = params

        # Unpack parameters and check bounds
        if not (0.0 <= m <= 0.20 and 0.1 <= p <= 0.8 and 0.05 <= t <= 0.35):
            self.history.append([self.eval_count, m, p, t, "Bounds", "Bounds", self.FAILURE_PENALTY])
            return self.FAILURE_PENALTY

        # Generate airfoil
        X, Y, _ = naca4_airfoil(m, p, t)
        
        # Check geometric box constraint
        if self.max_height is not None:
            physical_height = (np.max(Y) - np.min(Y)) * self.chord
            if physical_height > self.max_height:
                height_penalty = ((physical_height - self.max_height) * 1000) ** 2
                score = self.FAILURE_PENALTY + height_penalty
                print(f"| {self.eval_count:4d} | {m:.4f} | {p:.4f} | {t:.4f} | {'Box Fail':^9} | {'Box Fail':^9} | {score:.4e} |")
                self.history.append([self.eval_count, m, p, t, "Box Fail", "Box Fail", score])
                return score
        
        # Run XFOIL
        analysis = XFoilAnalysis(
            airfoil_name=f"opt_naca_{self.eval_count}",
            alpha=self.alpha,
            reynolds=self.reynolds,
            mach=self.mach,
            ncrit=self.ncrit
        )
        cl, cd, achieved_alpha = analysis.run_analysis(X, Y)
        
        # Calculate score
        if cl is not None and cd is not None and cd > 0:
            cl_error = ((cl - self.target_cl) * self.CL_ERROR_WEIGHT) ** 2
            cd_penalty = (max(0, cd - self.max_cd) * self.CD_PENALTY_WEIGHT) ** 2
            
            alpha_penalty = 0
            alpha_mismatch_flag = ""
            if abs(achieved_alpha - self.alpha) > 0.1:
                alpha_penalty = (abs(achieved_alpha - self.alpha) * self.ALPHA_PENALTY_WEIGHT) ** 2
                alpha_mismatch_flag = "*"

            score = cl_error + cd_penalty + alpha_penalty
            
            print(f"| {self.eval_count:4d} | {m:.4f} | {p:.4f} | {t:.4f} | {cl:.4f}{alpha_mismatch_flag} | {cd:.5f}{alpha_mismatch_flag} | {score:.4e} |")
            self.history.append([self.eval_count, m, p, t, cl, cd, score])
            return score
        else:
            # Penalize non-converging solutions
            # Give a slight preference to more reasonable shapes even if they fail
            thickness_penalty = ((0.12 - t)**2 * 1e5)
            camber_penalty = ((0.05 - m)**2 * 1e5)
            emergency_score = self.FAILURE_PENALTY + thickness_penalty + camber_penalty
            
            print(f"| {self.eval_count:4d} | {m:.4f} | {p:.4f} | {t:.4f} | {'Failed':^9} | {'Failed':^9} | {emergency_score:.4e} |")
            self.history.append([self.eval_count, m, p, t, "Failed", "Failed", emergency_score])
            return emergency_score

    def optimize(self, initial_guess, bounds, max_iter=50):
        """
        Runs the optimization process.

        Args:
            initial_guess (list): Initial values for [m, p, t].
            bounds (list of tuples): Bounds for [(m), (p), (t)].
            max_iter (int): Maximum number of iterations for the optimizer.

        Returns:
            scipy.optimize.OptimizeResult: The result object from scipy.minimize.
        """
        print("[+] Starting XFOIL-based Airfoil Optimization...")
        print("-" * 70)
        print(f"| {'Eval':^4} | {'m':^6} | {'p':^6} | {'t':^6} | {'Cl':^9} | {'Cd':^9} | {'Score':^10} |")
        print("-" * 70)

        self.eval_count = 0
        self.history = []

        result = minimize(
            self._objective_function,
            initial_guess,
            method='SLSQP',
            bounds=bounds,
            options={'disp': True, 'maxiter': max_iter, 'ftol': 1e-4, 'eps': 1e-4}
        )
        
        print("-" * 70)
        return result

    def get_optimization_history(self):
        return self.history
