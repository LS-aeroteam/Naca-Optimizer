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
    """Calculates the source and vortex influence coefficients."""
    num_panels = len(XC)
    I = np.zeros((num_panels, num_panels))
    J = np.zeros((num_panels, num_panels))
    K = np.zeros((num_panels, num_panels))
    L = np.zeros((num_panels, num_panels))
    
    for i in range(num_panels):
        for j in range(num_panels):
            if i == j:
                continue

            A = -(XC[i]-XB[j])*_cosd(PSI[j]) - (YC[i]-YB[j])*_sind(PSI[j])
            B = (XC[i]-XB[j])**2 + (YC[i]-YB[j])**2
            E = math.sqrt(max(B - A**2, 0)) # Ensure non-negative for sqrt

            # Vortex-induced normal and tangential velocities
            Cn_v = -_cosd(PSI[i]-PSI[j])
            Dn_v = (XC[i]-XB[j])*_cosd(PSI[i]) + (YC[i]-YB[j])*_sind(PSI[i])
            Ct_v = _sind(PSI[j]-PSI[i])
            Dt_v = (XC[i]-XB[j])*_sind(PSI[i]) - (YC[i]-YB[j])*_cosd(PSI[i])
            
            term1_v = 0.5 * Cn_v * math.log((S[j]**2 + 2*A*S[j] + B) / B)
            term2_v = ((Dn_v - A*Cn_v) / E) * (math.atan2(S[j]+A, E) - math.atan2(A, E)) if E != 0 else 0
            K[i, j] = term1_v + term2_v
            
            term1_l = 0.5*Ct_v*math.log((S[j]**2 + 2*A*S[j] + B)/B)
            term2_l = ((Dt_v - A*Ct_v) / E) * (math.atan2(S[j]+A, E) - math.atan2(A, E)) if E != 0 else 0
            L[i,j] = term1_l + term2_l

            # Source-induced normal and tangential velocities
            Cn_s = _sind(PSI[i]-PSI[j])
            Dn_s = -(XC[i]-XB[j])*_sind(PSI[i]) + (YC[i]-YB[j])*_cosd(PSI[i])
            
            term1_i = 0.5 * Cn_s * math.log((S[j]**2 + 2*A*S[j] + B) / B)
            term2_i = ((Dn_s - A*Cn_s) / E) * (math.atan2(S[j]+A, E) - math.atan2(A, E)) if E != 0 else 0
            I[i, j] = term1_i + term2_i
            
            Ct_s = -_cosd(PSI[i]-PSI[j])
            Dt_s = (XC[i]-XB[j])*_cosd(PSI[i]) + (YC[i]-YB[j])*_sind(PSI[i])

            term1_j = 0.5 * Ct_s * math.log((S[j]**2 + 2*A*S[j] + B) / B)
            term2_j = ((Dt_s - A*Ct_s) / E) * (math.atan2(S[j]+A, E) - math.atan2(A, E)) if E != 0 else 0
            J[i, j] = term1_j + term2_j

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

def run_panel_analysis(XB, YB, alpha_deg):
    """
    Runs a complete panel method analysis for a given airfoil geometry and angle of attack.

    Args:
        XB (np.ndarray): Airfoil X-coordinates (expected to be TE -> LE -> TE).
        YB (np.ndarray): Airfoil Y-coordinates.
        alpha_deg (float): Angle of attack in degrees.

    Returns:
        dict: A dictionary containing the results:
              'cl_potential' (float), 'Cp' (np.ndarray), 'XC' (np.ndarray), 'YC' (np.ndarray).
    """
    # Enforce clockwise order for the panel method calculations
    # naca4_airfoil returns counter-clockwise (TE -> Upper -> LE -> Lower -> TE)
    # We reverse it to be TE -> Lower -> LE -> Upper -> TE (Clockwise)
    XB = XB[::-1]
    YB = YB[::-1]
    
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
    
    print(f"    Potential Lift Coefficient (CL): {cl_potential:.4f}")
    
    return {
        'cl_potential': cl_potential,
        'Cp': Cp,
        'XC': XC,
        'YC': YC,
        'num_panels': len(XB) - 1
    }
