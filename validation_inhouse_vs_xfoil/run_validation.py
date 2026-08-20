import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "xfoil_viscous_optimizer")))

from xfoil_core.airfoil import naca4_airfoil
from xfoil_core.panel_method import run_panel_analysis
from xfoil_core.xfoil import XFoilAnalysis
from xfoil_core.pre_run_checks import perform_all_checks

def plot_validation_results(results_dict, alphas, results_subfolder):
    """Generates comparative Cl vs Alpha plots."""
    plt.figure(figsize=(10, 6))
    
    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    
    for idx, (profile_name, data) in enumerate(results_dict.items()):
        c = colors[idx % len(colors)]
        
        # Filter valid values
        valid_alphas = [d['Alpha'] for d in data if d['Cl_Xfoil'] is not None]
        cl_inhouse = [d['Cl_InHouse'] for d in data if d['Cl_Xfoil'] is not None]
        cl_xfoil = [d['Cl_Xfoil'] for d in data if d['Cl_Xfoil'] is not None]
        
        if not valid_alphas:
            continue
            
        plt.plot(valid_alphas, cl_inhouse, marker='o', linestyle='--', color=c, label=f'{profile_name} (In-House)')
        plt.plot(valid_alphas, cl_xfoil, marker='s', linestyle='-', color=c, label=f'{profile_name} (XFOIL)')

    plt.title('Lift Coefficient (Cl) vs Angle of Attack (Alpha) Comparison')
    plt.xlabel('Alpha (degrees)')
    plt.ylabel('Cl')
    plt.grid(True)
    plt.legend()
    
    os.makedirs(results_subfolder, exist_ok=True)
    plot_filename = os.path.join(results_subfolder, "validation_plot_cl.svg")
    plt.savefig(plot_filename)
    print(f"[i] Plot saved to: {plot_filename}")

def get_fluid_selection():
    fluids = {
        '1': {'name': 'Air (Standard SL)', 'viscosity': 1.46e-5, 'speed_of_sound': 340.3},
        '2': {'name': 'Water (20 degrees)', 'viscosity': 1.00e-6, 'speed_of_sound': 1482.0}
    }
    while True:
        print("\nSelect the operating fluid:")
        print("1. Air (Kinematic viscosity: 1.46e-5 m^2/s, Speed of sound: 340.3 m/s)")
        print("2. Water (Kinematic viscosity: 1.00e-6 m^2/s, Speed of sound: 1482.0 m/s)")
        choice = input("Choice (1 or 2) [default: 1]: ").strip()
        if not choice:
            return fluids['1']
        if choice in fluids:
            return fluids[choice]
        print("Error: invalid selection.")

def get_user_input():
    print("\n======================================================================")
    print("                SETUP FOR VALIDATION")
    print("======================================================================")
    
    fluid = get_fluid_selection()
    
    try:
        alpha_min = float(input("\nEnter minimum alpha (degrees) [default: -4]: ") or -4)
        alpha_max = float(input("Enter maximum alpha (degrees) [default: 10]: ") or 10)
        alpha_step = float(input("Enter alpha step (degrees) [default: 2]: ") or 2)
        alphas = np.arange(alpha_min, alpha_max + (alpha_step/2), alpha_step)
        
        speed = float(input("Enter design speed (m/s) [e.g., 50]: ") or 50)
        chord = float(input("Enter airfoil chord (m) [e.g., 1.0]: ") or 1.0)
        
        reynolds = (speed * chord) / fluid['viscosity']
        mach = speed / fluid['speed_of_sound']
        print(f"\n [i] Calculated Conditions:")
        print(f"    - Reynolds Number: {reynolds:,.0f}")
        print(f"    - Mach Number: {mach:.3f}")

        num_points = int(input("Enter number of points/panels [default: 160]: ") or 160)
        return alphas, reynolds, mach, num_points
    except ValueError:
        print("\n[!] Invalid input. Using default values.")
        return np.arange(-4, 12, 2), 1e6, 0.0, 160

def main():
    print("=========================================================")
    print("       IN-HOUSE SOLVER vs XFOIL VALIDATION SCRIPT")
    print("=========================================================")
    
    # Run preliminary checks (ensures XFOIL is present)
    perform_all_checks()

    alphas, reynolds, mach, num_points = get_user_input()

    # TEST MATRIX
    # We test some common NACA profiles: symmetric, mildly cambered, and highly cambered
    profiles = [
        (0.0, 0.0, 0.12, "0012"),  # Symmetric
        (0.02, 0.4, 0.12, "2412"), # Standard cambered
        (0.04, 0.4, 0.12, "4412")  # Highly cambered
    ]
    
    # Lists and dictionaries to collect data
    all_results = []
    plot_data = {}

    print(f"Testing {len(profiles)} profiles and {len(alphas)} angles of attack (Total runs: {len(profiles)*len(alphas)})")
    
    for m, p, t, name in profiles:
        print(f"\n--- Processing Profile NACA {name} ---")
        profile_name = f"NACA {name}"
        plot_data[profile_name] = []
        
        # 1. Geometry Generation
        X_panel, Y_panel, _ = naca4_airfoil(m, p, t, num_points=int(num_points/2)+1)
        
        # 2. XFOIL Analyzer Initialization
        xfoil_analyzer = XFoilAnalysis(airfoil_name=f"naca_{name}", alpha=0.0, reynolds=reynolds, mach=mach)
        
        for alpha in alphas:
            print(f"  > Alpha = {alpha:4.1f} deg ...", end=" ", flush=True)
            
            # --- In-House Solver (Panel Method) ---
            # run_panel_analysis in the Xfoil wrapper does not support verbose=False, will print to screen
            try:
                panel_res = run_panel_analysis(X_panel, Y_panel, float(alpha))
                cl_inhouse = panel_res['cl_potential']
            except Exception as e:
                cl_inhouse = 0.0
                print(f"[ERR Panels: {e}]", end=" ")
            
            # --- XFOIL ---
            xfoil_analyzer.alpha = float(alpha)
            cl_xfoil, cd_xfoil, achieved_alpha = xfoil_analyzer.run_analysis(X_panel, Y_panel)
            
            # --- Error Calculation ---
            if cl_xfoil is not None and cl_xfoil != 0:
                err_cl = abs((cl_inhouse - cl_xfoil) / cl_xfoil) * 100
            else:
                err_cl = None
                
            # Print quick result
            status = "OK" if cl_xfoil is not None else "XFOIL TIMEOUT"
            cl_xfoil_str = f"{cl_xfoil:.4f}" if cl_xfoil is not None else "N/A"
            err_str = f"{err_cl:.1f}%" if err_cl is not None else "N/A"
            
            print(f"[{status}] In-House Cl: {cl_inhouse:.4f} | XFOIL Cl: {cl_xfoil_str} | Err: {err_str}")
            
            res_dict = {
                'Profile': profile_name,
                'Alpha': alpha,
                'Cl_InHouse': cl_inhouse,
                'Cl_Xfoil': cl_xfoil,
                'Error_Cl_%': err_cl
            }
            
            all_results.append(res_dict)
            plot_data[profile_name].append(res_dict)

    # --- Save CSV ---
    base_dir = os.path.dirname(os.path.abspath(__file__))
    results_subfolder_name = f"Validation_Re{int(reynolds)}_Mach{mach:.3f}"
    results_subfolder = os.path.join(base_dir, "Results", results_subfolder_name)
    os.makedirs(results_subfolder, exist_ok=True)
    csv_file = os.path.join(results_subfolder, "validation_results.csv")
    print(f"\n[+] Saving test results to '{csv_file}'...")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Profile', 'Alpha', 'Cl_InHouse', 'Cl_Xfoil', 'Error_Cl_%'])
        writer.writeheader()
        writer.writerows(all_results)
        
    # --- Plot Generation ---
    print(f"[+] Generating comparative plots...")
    plot_validation_results(plot_data, alphas, results_subfolder)
        
    print("\n[+] Validation completed successfully!")

if __name__ == "__main__":
    main()
