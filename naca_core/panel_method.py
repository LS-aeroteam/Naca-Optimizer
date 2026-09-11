import numpy as np
import math

def _cosd(deg): return np.cos(np.radians(deg))
def _sind(deg): return np.sin(np.radians(deg))
def _atan2d(y, x): 
    res = math.degrees(math.atan2(y, x))
    return res if res >= 0 else res + 360

def _calculate_geometry_parameters(XB, YB):
    """Calculates panel control points, lengths, and angles."""
    num_panels = len(XB) - 1
    XC = np.zeros(num_panels)
    YC = np.zeros(num_panels)
    S = np.zeros(num_panels)
    PSI = np.zeros(num_panels)
    
    for i in range(num_panels):
        XC[i] = (XB[i] + XB[i+1]) / 2
        YC[i] = (YB[i] + YB[i+1]) / 2
        S[i] = math.hypot(XB[i+1] - XB[i], YB[i+1] - YB[i])
        PSI[i] = _atan2d(YB[i+1] - YB[i], XB[i+1] - XB[i])
        
    return XC, YC, S, PSI

def _calculate_influence_coefficients(XC, YC, XB, YB, PSI, S):
    """
    Calculates the source and vortex influence coefficients.

    Row i = control point where the velocity is evaluated, column j = panel that induces it.
    The formulas are the same as the original double loop, evaluated for all (i, j) pairs at
    once with NumPy arrays. The diagonal (i == j) is zero, as in the loop.
    """
    num_panels = len(XC)
    XC = np.asarray(XC, dtype=float)[:, None]
    YC = np.asarray(YC, dtype=float)[:, None]
    XBj = np.asarray(XB, dtype=float)[None, :num_panels]
    YBj = np.asarray(YB, dtype=float)[None, :num_panels]
    PSI = np.asarray(PSI, dtype=float)
    PSIi, PSIj = PSI[:, None], PSI[None, :]
    Sj = np.asarray(S, dtype=float)[None, :]

    dx = XC - XBj
    dy = YC - YBj
    cos_i, sin_i = _cosd(PSIi), _sind(PSIi)
    cos_j, sin_j = _cosd(PSIj), _sind(PSIj)

    A = -dx * cos_j - dy * sin_j
    B = dx**2 + dy**2
    E = np.sqrt(np.maximum(B - A**2, 0))  # Ensure non-negative for sqrt

    off_diag = ~np.eye(num_panels, dtype=bool)
    log_term = np.zeros_like(B)
    log_term[off_diag] = np.log(((Sj**2 + 2*A*Sj + B) / B)[off_diag])
    atan_term = np.arctan2(Sj + A, E) - np.arctan2(A, E)
    # Where E == 0 the second term is zero (as in the loop): divide only where E != 0
    atan_over_E = np.divide(atan_term, E, out=np.zeros_like(E), where=(E != 0) & off_diag)

    # Vortex-induced normal and tangential velocities
    Cn_v = -_cosd(PSIi - PSIj)
    Dn_v = dx * cos_i + dy * sin_i
    Ct_v = _sind(PSIj - PSIi)
    Dt_v = dx * sin_i - dy * cos_i
    K = 0.5 * Cn_v * log_term + (Dn_v - A * Cn_v) * atan_over_E
    L = 0.5 * Ct_v * log_term + (Dt_v - A * Ct_v) * atan_over_E

    # Source-induced normal and tangential velocities
    Cn_s = _sind(PSIi - PSIj)
    Dn_s = -dx * sin_i + dy * cos_i
    Ct_s = -_cosd(PSIi - PSIj)
    Dt_s = dx * cos_i + dy * sin_i
    I = 0.5 * Cn_s * log_term + (Dn_s - A * Cn_s) * atan_over_E
    J = 0.5 * Ct_s * log_term + (Dt_s - A * Ct_s) * atan_over_E

    for M in (I, J, K, L):
        np.fill_diagonal(M, 0.0)

    return I, J, K, L

def _solve_linear_system(I, J, K, L, PSI, alpha_deg):
    """Builds and solves the linear system for source and vortex strengths."""
    num_panels = len(PSI)
    A_mat = np.zeros((num_panels + 1, num_panels + 1))
    
    A_mat[:num_panels, :num_panels] = I
    np.fill_diagonal(A_mat, np.pi)
    
    A_mat[:num_panels, num_panels] = np.sum(K, axis=1)
    A_mat[num_panels, :num_panels] = J[0, :] + J[num_panels-1, :]
    A_mat[num_panels, num_panels] = np.sum(L[0, :] + L[num_panels-1, :]) - 2 * np.pi
    
    b = np.zeros(num_panels + 1)
    b[:num_panels] = -2 * np.pi * _sind(alpha_deg - PSI)
    b[num_panels] = -2 * np.pi * (_cosd(alpha_deg - PSI[0]) + _cosd(alpha_deg - PSI[num_panels-1]))
    
    res = np.linalg.solve(A_mat, b)
    lambda_src = res[:-1]
    gamma = res[-1]
    
    return lambda_src, gamma

def _calculate_surface_velocities(lambda_src, gamma, J, L, PSI, alpha_deg):
    """Calculates tangential velocities and pressure coefficients."""
    num_panels = len(PSI)
    Vt = np.zeros(num_panels)
    
    for i in range(num_panels):
        term1 = _cosd(alpha_deg - PSI[i])
        term2 = (1 / (2 * np.pi)) * np.sum(lambda_src * J[i, :])
        term3 = -gamma / 2
        term4 = (gamma / (2 * np.pi)) * np.sum(L[i, :])
        Vt[i] = term1 + term2 + term3 + term4
        
    Cp = 1 - (Vt**2)
    return Vt, Cp

def surface_distributions(panel_results):
    """
    Splits Cp into upper and lower surface and aligns them on the same x/c.

    For cambered airfoils the upper and lower control points with the same index
    are not at the same x/c, so the lower-surface Cp is linearly interpolated at
    the upper-surface control points before computing Delta Cp.

    Args:
        panel_results (dict): The results dictionary from `run_panel_analysis`.

    Returns:
        tuple: (x_upper, cp_upper, cp_lower, delta_cp), ordered from LE to TE,
            where cp_lower is evaluated at x_upper and delta_cp = cp_lower - cp_upper.
    """
    n_half = panel_results['num_panels'] // 2
    XC = np.asarray(panel_results['XC'])
    Cp = np.asarray(panel_results['Cp'])

    x_upper, cp_upper = XC[n_half:], Cp[n_half:]
    x_lower, cp_lower = XC[:n_half], Cp[:n_half]

    order = np.argsort(x_lower)  # np.interp needs increasing x
    cp_lower_at_upper = np.interp(x_upper, x_lower[order], cp_lower[order])

    return x_upper, cp_upper, cp_lower_at_upper, cp_lower_at_upper - cp_upper


def run_panel_analysis(XB, YB, alpha_deg, verbose=True):
    """
    Runs a complete panel method analysis for a given airfoil geometry and angle of attack.

    Args:
        XB (np.ndarray): Airfoil X-coordinates (expected to be TE -> LE -> TE).
        YB (np.ndarray): Airfoil Y-coordinates.
        alpha_deg (float): Angle of attack in degrees.
        verbose (bool): Whether to print output.

    Returns:
        dict: A dictionary containing the results:
              'cl_potential' (float), 'Cp' (np.ndarray), 'XC' (np.ndarray), 'YC' (np.ndarray).
    """
    # Enforce clockwise order for the panel method calculations
    # naca4_airfoil returns counter-clockwise (TE -> Upper -> LE -> Lower -> TE)
    # We reverse it to be TE -> Lower -> LE -> Upper -> TE (Clockwise)
    XB = XB[::-1]
    YB = YB[::-1]
    
    if verbose:
        print("\n[+] Starting Panel Method Simulation (Sources + Vortices)...")
    
    # 1. Calculate geometry parameters
    XC, YC, S, PSI = _calculate_geometry_parameters(XB, YB)
    
    # 2. Calculate influence coefficients
    I, J, K, L = _calculate_influence_coefficients(XC, YC, XB, YB, PSI, S)
    
    # 3. Solve for source and vortex strengths
    lambda_src, gamma = _solve_linear_system(I, J, K, L, PSI, alpha_deg)
    
    # 4. Calculate surface velocities and Cp
    Vt, Cp = _calculate_surface_velocities(lambda_src, gamma, J, L, PSI, alpha_deg)
    
    # 5. Calculate total lift coefficient
    perimeter = np.sum(S)
    gamma_total = gamma * perimeter
    cl_potential = -2 * gamma_total  # Based on Kutta-Joukowski theorem L = -rho*V*Gamma
    
    if verbose:
        print(f"    Potential Lift Coefficient (CL): {cl_potential:.4f}")
    
    return {
        'cl_potential': cl_potential,
        'Cp': Cp,
        'XC': XC,
        'YC': YC,
        'num_panels': len(XB) - 1
    }
