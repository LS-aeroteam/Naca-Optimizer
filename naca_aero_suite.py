import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import subprocess
import os
import csv
import math
from scipy.optimize import minimize

# ==============================================================================
# 1. FUNZIONI DI GENERAZIONE PROFILO E UTILITY (XFOIL)
# ==============================================================================

def naca4(m_param, p_param, t_param, c=1.0, n=100):
    """Genera le coordinate del profilo NACA 4 cifre con spaziatura cosinusoidale per XFOIL."""
    beta = np.linspace(0, np.pi, n)
    x = c * (0.5 * (1 - np.cos(beta)))
    yt = 5 * t_param * c * (0.2969 * np.sqrt(x/c) - 0.1260 * (x/c) - 0.3516 * (x/c)**2 + 0.2843 * (x/c)**3 - 0.1015 * (x/c)**4)

    if p_param == 0 or m_param == 0:
        xu, yu = x, yt
        xl, yl = x, -yt
    else:
        yc = np.zeros_like(x)
        dyc_dx = np.zeros_like(x)
        
        front_x = x[x < p_param * c]
        back_x = x[x >= p_param * c]

        if len(front_x) > 0:
            yc_front = (m_param / p_param**2) * (2 * p_param * (front_x / c) - (front_x / c)**2)
            dyc_dx_front = (2 * m_param / p_param**2) * (p_param - front_x / c)
            yc[:len(front_x)] = yc_front
            dyc_dx[:len(front_x)] = dyc_dx_front

        if len(back_x) > 0:
            yc_back = (m_param / (1 - p_param)**2) * ((1 - 2 * p_param) + 2 * p_param * (back_x / c) - (back_x / c)**2)
            dyc_dx_back = (2 * m_param / (1 - p_param)**2) * (p_param - back_x / c)
            yc[len(front_x):] = yc_back
            dyc_dx[len(front_x):] = dyc_dx_back
            
        theta = np.arctan(dyc_dx)
        xu = x - yt * np.sin(theta)
        yu = yc + yt * np.cos(theta)
        xl = x + yt * np.sin(theta)
        yl = yc - yt * np.cos(theta)

    X = np.concatenate((np.flip(xu), xl[1:]))
    Y = np.concatenate((np.flip(yu), yl[1:]))
    return X, Y, (xu, yu, xl, yl)

def save_airfoil_to_file(X, Y, filename):
    with open(filename, "w") as f:
        for i in range(len(X)):
            f.write(f"{X[i]:.6f} {Y[i]:.6f}\n")

# ==============================================================================
# 2. WRAPPER CFD E OTTIMIZZAZIONE (XFOIL + SLSQP)
# ==============================================================================

def get_xfoil_executable():
    """Riconosce automaticamente l'ambiente (Windows nativo o Linux/WSL)."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if os.name == 'nt':
        return os.path.join(script_dir, "xfoil.exe")
    else:
        return "xfoil" # Assumiamo sia nel PATH o nella cartella in ambiente Linux/WSL

def run_xfoil_analysis(airfoil_file, alpha, Re, Mach=0.0, ncrit=9.0):
    xfoil_input_file = "xfoil_input.in"
    polar_file = "polar.dat"

    for f in [xfoil_input_file, polar_file]:
        if os.path.exists(f): os.remove(f)

    with open(xfoil_input_file, "w") as f:
        f.write(f"LOAD {airfoil_file}\n")
        f.write("PANE\nOPER\nVPAR\n")
        f.write(f"N {ncrit}\n\n")
        f.write("ALFA 0\n") 
        f.write(f"Visc {Re}\n")
        f.write(f"Mach {Mach}\n")
        f.write("ITER 500\nPACC\n")
        f.write(f"{polar_file}\n\n") 
        
        if alpha == 0.0:
            f.write("ALFA 0.0\n")
        else:
            step = 1.0 if alpha > 0 else -1.0
            f.write(f"ASEQ 0.0 {alpha - step/2} {step}\n")
            f.write(f"ALFA {alpha}\n")
        f.write("\nQUIT\n")

    xfoil_exe = get_xfoil_executable()

    try:
        with open(xfoil_input_file, "r") as stdin_file:
            subprocess.run([xfoil_exe], stdin=stdin_file, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=30)
    except Exception:
        pass

    cl, cd, achieved_alpha = None, None, None
    try:
        with open(polar_file, "r") as f:
            lines = [line for line in f if not line.startswith("#") and len(line.strip()) > 0]
            best_diff = 1e6
            for line in lines:
                data = line.split()
                if len(data) >= 3:
                    try:
                        a_val, cl_val, cd_val = float(data[0]), float(data[1]), float(data[2])
                        if abs(a_val - alpha) < best_diff:
                            best_diff = abs(a_val - alpha)
                            achieved_alpha, cl, cd = a_val, cl_val, cd_val
                    except ValueError:
                        pass
    except (IOError, IndexError):
        pass
    
    for f in [xfoil_input_file, polar_file, airfoil_file]:
        if os.path.exists(f): os.remove(f)
            
    return cl, cd, achieved_alpha

def objective_function(params, Re, alpha, target_cl, max_cd, mach, ncrit):
    objective_function.eval_count += 1
    eval_count = objective_function.eval_count
    m, p, t = params
    
    if not (0.0 <= m <= 0.20 and 0.1 <= p <= 0.8 and 0.05 <= t <= 0.35):
        objective_function.history.append([eval_count, m, p, t, "Bounds", "Bounds", 1e6])
        return 1e6

    airfoil_name = f"temp_naca_{eval_count}.dat"
    X, Y, _ = naca4(m, p, t)
    save_airfoil_to_file(X, Y, airfoil_name)
    cl, cd, achieved_alpha = run_xfoil_analysis(airfoil_name, alpha, Re, mach, ncrit)
    
    if cl is not None and cd is not None and cd > 0:
        cl_error = ((cl - target_cl) * 10) ** 2
        cd_penalty = (max(0, cd - max_cd) * 1000) ** 2
        
        if abs(achieved_alpha - alpha) > 0.1:
            alpha_penalty = (abs(achieved_alpha - alpha) * 1000) ** 2
            score = cl_error + cd_penalty + alpha_penalty
            print(f"| {eval_count:4d} | {m:.4f} | {p:.4f} | {t:.4f} | {cl:.4f}* | {cd:.5f}* | {score:.4e} |")
            objective_function.history.append([eval_count, m, p, t, cl, cd, score])
            return score
            
        score = cl_error + cd_penalty
        print(f"| {eval_count:4d} | {m:.4f} | {p:.4f} | {t:.4f} |  {cl:.4f}  |  {cd:.5f} | {score:.4e} |")
        objective_function.history.append([eval_count, m, p, t, cl, cd, score])
        return score
    else:
        emergency = 1e6 + ((0.12 - t)**2 * 1e5) + ((0.05 - m)**2 * 1e5)
        print(f"| {eval_count:4d} | {m:.4f} | {p:.4f} | {t:.4f} | Separated | Separated | {emergency:.4e} |")
        objective_function.history.append([eval_count, m, p, t, "Separated", "Separated", emergency])
        return emergency

# ==============================================================================
# 3. ANALISI AVANZATA (METODO DEI PANNELLI TRADOTTO DA MATLAB)
# ==============================================================================

def cosd(deg): return math.cos(math.radians(deg))
def sind(deg): return math.sin(math.radians(deg))
def atan2d(y, x): 
    res = math.degrees(math.atan2(y, x))
    return res if res >= 0 else res + 360

def run_panel_method(m, p, t, alpha_deg, nPanels=160, plot_streamlines=True, plot_cp=True):
    print("\n[+] Avvio simulazione Metodo dei Pannelli (Sorgenti + Vortici)...")
    
    # --- Generazione Geometria Rigorosa (MATLAB Logic) ---
    nHalf = int(nPanels / 2)
    beta = np.linspace(0, np.pi, nHalf + 1)
    x = (1 - np.cos(beta)) / 2
    
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    yc = np.zeros_like(x)
    dyc = np.zeros_like(x)
    
    if m != 0:
        for i, xi in enumerate(x):
            if xi < p:
                yc[i] = (m / p**2) * (2 * p * xi - xi**2)
                dyc[i] = (2 * m / p**2) * (p - xi)
            else:
                yc[i] = m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * xi - xi**2)
                dyc[i] = 2 * m / (1 - p)**2 * (p - xi)
                
    theta = np.arctan(dyc)
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    
    # Punti in senso orario: ventre (BU->BA) + dorso (BA->BU)
    XB = np.concatenate((xl[::-1], xu[1:]))
    YB = np.concatenate((yl[::-1], yu[1:]))
    nPanels = len(XB) - 1
    
    XC = np.zeros(nPanels)
    YC = np.zeros(nPanels)
    S = np.zeros(nPanels)
    PSI = np.zeros(nPanels)
    
    for i in range(nPanels):
        XC[i] = (XB[i] + XB[i+1]) / 2
        YC[i] = (YB[i] + YB[i+1]) / 2
        S[i] = math.hypot(XB[i+1] - XB[i], YB[i+1] - YB[i])
        PSI[i] = atan2d(YB[i+1] - YB[i], XB[i+1] - XB[i])
        
    # --- Matrici di Influenza ---
    K = np.zeros((nPanels, nPanels))
    L = np.zeros((nPanels, nPanels))
    I = np.zeros((nPanels, nPanels))
    J = np.zeros((nPanels, nPanels))
    
    for i in range(nPanels):
        for j in range(nPanels):
            if j != i:
                A = -(XC[i]-XB[j])*cosd(PSI[j]) - (YC[i]-YB[j])*sind(PSI[j])
                B = (XC[i]-XB[j])**2 + (YC[i]-YB[j])**2
                Cn = -cosd(PSI[i]-PSI[j])
                Dn = (XC[i]-XB[j])*cosd(PSI[i]) + (YC[i]-YB[j])*sind(PSI[i])
                Ct = sind(PSI[j]-PSI[i])
                Dt = (XC[i]-XB[j])*sind(PSI[i]) - (YC[i]-YB[j])*cosd(PSI[i])
                E = math.sqrt(max(B - A**2, 0))
                
                term1_v = 0.5*Cn*math.log((S[j]**2 + 2*A*S[j] + B)/B)
                term2_v = ((Dn-A*Cn)/E) * (math.atan2(S[j]+A, E) - math.atan2(A, E)) if E != 0 else 0
                K[i,j] = term1_v + term2_v
                
                term1_l = 0.5*Ct*math.log((S[j]**2 + 2*A*S[j] + B)/B)
                term2_l = ((Dt-A*Ct)/E) * (math.atan2(S[j]+A, E) - math.atan2(A, E)) if E != 0 else 0
                L[i,j] = term1_l + term2_l
                
                Cn_s = sind(PSI[i]-PSI[j])
                Dn_s = -(XC[i]-XB[j])*sind(PSI[i]) + (YC[i]-YB[j])*cosd(PSI[i])
                Ct_s = -cosd(PSI[i]-PSI[j])
                Dt_s = (XC[i]-XB[j])*cosd(PSI[i]) + (YC[i]-YB[j])*sind(PSI[i])
                
                term1_i = 0.5*Cn_s*math.log((S[j]**2 + 2*A*S[j] + B)/B)
                term2_i = ((Dn_s-A*Cn_s)/E) * (math.atan2(S[j]+A, E) - math.atan2(A, E)) if E!=0 else 0
                I[i,j] = term1_i + term2_i
                
                term1_j = 0.5*Ct_s*math.log((S[j]**2 + 2*A*S[j] + B)/B)
                term2_j = ((Dt_s-A*Ct_s)/E) * (math.atan2(S[j]+A, E) - math.atan2(A, E)) if E!=0 else 0
                J[i,j] = term1_j + term2_j

    # --- Sistema Lineare e Condizione di Kutta ---
    A_mat = np.zeros((nPanels+1, nPanels+1))
    for i in range(nPanels):
        for j in range(nPanels):
            if i == j:
                A_mat[i,j] = np.pi
            else:
                A_mat[i,j] = I[i,j]
                
    for i in range(nPanels):
        A_mat[i, nPanels] = np.sum(K[i,:])
        
    for j in range(nPanels):
        A_mat[nPanels, j] = J[0, j] + J[nPanels-1, j]
        
    A_mat[nPanels, nPanels] = np.sum(L[0,:] + L[nPanels-1,:]) - 2*np.pi
    
    b = np.zeros(nPanels+1)
    for i in range(nPanels):
        b[i] = -2*np.pi * sind(alpha_deg - PSI[i])
    b[nPanels] = -2*np.pi * (cosd(alpha_deg - PSI[0]) + cosd(alpha_deg - PSI[nPanels-1]))
    
    resArr = np.linalg.solve(A_mat, b)
    lambda_src = resArr[:-1]
    gamma = resArr[-1]
    
    Perimeter = np.sum(S)
    GAMMA_total = gamma * Perimeter
    CL_inviscid = -2 * GAMMA_total
    
    print(f"CL Calcolato (Inviscido/Teorico): {CL_inviscid:.4f}")
    
    # --- Calcolo Velocità Tangenziali e Cp ---
    Vt = np.zeros(nPanels)
    Cp = np.zeros(nPanels)
    for i in range(nPanels):
        term1 = cosd(alpha_deg - PSI[i])
        term2 = (1/(2*np.pi)) * np.sum(lambda_src * J[i,:])
        term3 = -gamma/2
        term4 = (gamma/(2*np.pi)) * np.sum(L[i,:])
        Vt[i] = term1 + term2 + term3 + term4
        Cp[i] = 1 - (Vt[i]**2)

    # --- Plot Cp ---
    if plot_cp:
        xLower = XC[:nHalf]
        CpLower = Cp[:nHalf]
        xUpper = XC[nHalf:]
        CpUpper = Cp[nHalf:]
        
        plt.figure(figsize=(10, 5))
        plt.plot(xUpper, CpUpper, 'b-', linewidth=1.5, label='Dorso (Upper)')
        plt.plot(xLower, CpLower, 'r-', linewidth=1.5, label='Ventre (Lower)')
        plt.gca().invert_yaxis()
        plt.xlabel('x/c')
        plt.ylabel('$C_p$')
        plt.title(f'Distribuzione $C_p$ Metodo Pannelli - $\\alpha$ = {alpha_deg}°')
        plt.grid(True)
        plt.legend()
        plt.show()

    # --- Plot Streamlines ---
    if plot_streamlines:
        print("[+] Calcolo campo di moto per le streamlines. Attendi...")
        nGridX, nGridY = 100, 100
        xVals = [np.min(XB)-0.5, np.max(XB)+0.5]
        yVals = [np.min(YB)-0.3, np.max(YB)+0.3]
        
        Xgrid = np.linspace(xVals[0], xVals[1], nGridX)
        Ygrid = np.linspace(yVals[0], yVals[1], nGridY)
        XX, YY = np.meshgrid(Xgrid, Ygrid)
        
        Vx = np.zeros((nGridY, nGridX))
        Vy = np.zeros((nGridY, nGridX))
        
        # Uso MPath per velocizzare "inpolygon"
        airfoil_path = mpath.Path(np.column_stack((XB, YB)))
        points = np.column_stack((XX.flatten(), YY.flatten()))
        inside_mask = airfoil_path.contains_points(points).reshape(nGridY, nGridX)
        
        for n in range(nGridY):
            for m in range(nGridX):
                if inside_mask[n, m]:
                    Vx[n, m] = 0
                    Vy[n, m] = 0
                    continue
                    
                XP, YP = XX[n, m], YY[n, m]
                Mx, My, Nx, Ny = np.zeros(nPanels), np.zeros(nPanels), np.zeros(nPanels), np.zeros(nPanels)
                
                for j in range(nPanels):
                    A = -(XP-XB[j])*cosd(PSI[j]) - (YP-YB[j])*sind(PSI[j])
                    B = (XP-XB[j])**2 + (YP-YB[j])**2
                    E = math.sqrt(max(0, B - A**2))
                    
                    Cx, Dx = -cosd(PSI[j]), XP - XB[j]
                    Cy, Dy = -sind(PSI[j]), YP - YB[j]
                    
                    if E != 0:
                        term1_mx = 0.5*Cx*math.log((S[j]**2+2*A*S[j]+B)/B)
                        term2_mx = ((Dx-A*Cx)/E)*(math.atan2(S[j]+A, E) - math.atan2(A, E))
                        Mx[j] = term1_mx + term2_mx
                        
                        term1_my = 0.5*Cy*math.log((S[j]**2+2*A*S[j]+B)/B)
                        term2_my = ((Dy-A*Cy)/E)*(math.atan2(S[j]+A, E) - math.atan2(A, E))
                        My[j] = term1_my + term2_my
                        
                        Cx_n, Dx_n = sind(PSI[j]), -(YP-YB[j])
                        Cy_n, Dy_n = -cosd(PSI[j]), XP-XB[j]
                        
                        term1_nx = 0.5*Cx_n*math.log((S[j]**2+2*A*S[j]+B)/B)
                        term2_nx = ((Dx_n-A*Cx_n)/E)*(math.atan2(S[j]+A, E) - math.atan2(A, E))
                        Nx[j] = term1_nx + term2_nx
                        
                        term1_ny = 0.5*Cy_n*math.log((S[j]**2+2*A*S[j]+B)/B)
                        term2_ny = ((Dy_n-A*Cy_n)/E)*(math.atan2(S[j]+A, E) - math.atan2(A, E))
                        Ny[j] = term1_ny + term2_ny

                Vx[n, m] = cosd(alpha_deg) + np.sum(lambda_src*Mx/(2*np.pi)) + np.sum(gamma*Nx/(2*np.pi))
                Vy[n, m] = sind(alpha_deg) + np.sum(lambda_src*My/(2*np.pi)) + np.sum(gamma*Ny/(2*np.pi))
                
        plt.figure(figsize=(12, 6))
        plt.fill(XB, YB, 'w', edgecolor='k', linewidth=1.5, zorder=5)
        
        # Genero streamlines dense
        start_y = np.linspace(yVals[0], yVals[1], 25)
        start_x = np.full_like(start_y, xVals[0])
        seed_points = np.column_stack([start_x, start_y])
        
        plt.streamplot(XX, YY, Vx, Vy, color=np.sqrt(Vx**2 + Vy**2), cmap='viridis', 
                       linewidth=1, density=2, start_points=seed_points, zorder=1)
        
        plt.title('Streamlines e Campo di Velocità', fontsize=16, fontweight='bold')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.grid(True)
        plt.xlim(xVals)
        plt.ylim(yVals)
        plt.show()

# ==============================================================================
# MAIN EXECUTION BLOCK
# ==============================================================================

if __name__ == "__main__":
    print("======================================================================")
    print(" AERODYNAMIC SUITE: INVERSE DESIGN & INVISCID ANALYSIS")
    print("======================================================================")
    
    fluids = {
        '1': {'name': 'Aria (Standard SL)', 'vis': 1.46e-5, 'sos': 340.3},
        '2': {'name': 'Acqua (20 gradi)', 'vis': 1.00e-6, 'sos': 1482.0}
    }
    
    fluid = fluids['1'] # Impostato fisso per comodità, espandibile come prima
    speed = float(input("Inserisci velocità di progetto (m/s) [es. 20]: ") or 20)
    chord = float(input("Inserisci corda (m) [es. 1.0]: ") or 1.0)
    TARGET_REYNOLDS = (speed * chord) / fluid['vis']
    TARGET_MACH = speed / fluid['sos']
    
    TARGET_ALPHA = float(input("\nInserisci angolo di attacco (gradi) [es. 3]: ") or 3.0)
    TARGET_CL = float(input("Inserisci Cl target [es. 0.6]: ") or 0.6)
    MAX_CD = float(input("Inserisci Cd massimo tollerato [es. 0.01]: ") or 0.01)
    
    initial_guess = [0.01, 0.4, 0.12]
    bounds = [(0.0, 0.20), (0.1, 0.8), (0.05, 0.35)]

    print("\n[+] Avvio Ottimizzazione XFOIL...")
    objective_function.eval_count = 0
    objective_function.history = []

    result = minimize(
        objective_function, initial_guess,
        args=(TARGET_REYNOLDS, TARGET_ALPHA, TARGET_CL, MAX_CD, TARGET_MACH, 9.0),
        method='SLSQP', bounds=bounds,
        options={'disp': True, 'maxiter': 50, 'ftol': 1e-4, 'eps': 1e-4}
    )

    if result.success or result.nfev > 0:
        opt_m, opt_p, opt_t = result.x
        naca_opt_str = f"{int(round(opt_m*100))}{int(round(opt_p*10))}{int(round(opt_t*100)):02d}"
        print(f"\n[+] PROFILO OTTIMIZZATO TROVATO: NACA {naca_opt_str}")
        
        # --- Modulo Metodo dei Pannelli ---
        ans = input("\nVuoi eseguire l'analisi avanzata (Metodo Pannelli) sul profilo ottimizzato? (s/n): ")
        if ans.lower().strip() == 's':
            n_pan = int(input("Inserisci il numero di pannelli [default: 160]: ") or 160)
            run_panel_method(opt_m, opt_p, opt_t, TARGET_ALPHA, nPanels=n_pan)
    else:
        print("Ottimizzazione fallita.")