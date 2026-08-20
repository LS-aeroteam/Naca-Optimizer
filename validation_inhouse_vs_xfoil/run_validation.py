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

def plot_validation_results(results_dict, alphas):
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
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    plot_filename = os.path.join(base_dir, "validation_plot_cl.svg")
    plt.savefig(plot_filename)
    print(f"[i] Plot saved to: {plot_filename}")

def main():
    print("=========================================================")
    print("       IN-HOUSE SOLVER vs XFOIL VALIDATION SCRIPT")
    print("=========================================================")
    
    # Run preliminary checks (ensures XFOIL is present)
    perform_all_checks()

    # TEST MATRIX
    # We test some common NACA profiles: symmetric, mildly cambered, and highly cambered
    profiles = [
        (0.0, 0.0, 0.12, "0012"),  # Symmetric
        (0.02, 0.4, 0.12, "2412"), # Standard cambered
        (0.04, 0.4, 0.12, "4412")  # Highly cambered
    ]
    
    # Range of angles of attack
    alphas = np.arange(-4, 12, 2)  # -4, -2, 0, 2, 4, 6, 8, 10
    
    # Flow conditions for XFOIL
    reynolds = 1e6
    mach = 0.0
    num_points = 160 # Number of points/panels
    
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
            print(f"  > Alpha = {alpha:2d} deg ...", end=" ", flush=True)
            
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
    csv_file = os.path.join(base_dir, "validation_results.csv")
    print(f"\n[+] Saving test results to '{csv_file}'...")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Profile', 'Alpha', 'Cl_InHouse', 'Cl_Xfoil', 'Error_Cl_%'])
        writer.writeheader()
        writer.writerows(all_results)
        
    # --- Plot Generation ---
    print(f"[+] Generating comparative plots...")
    plot_validation_results(plot_data, alphas)
        
    print("\n[+] Validation completed successfully!")

if __name__ == "__main__":
    main()
