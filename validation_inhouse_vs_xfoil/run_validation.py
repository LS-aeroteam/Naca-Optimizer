import os
import sys
import csv

# Make the shared 'naca_core' package (repository root) importable from this folder
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Only the standard library is imported here: numpy/matplotlib and the solver modules
# are imported inside the functions, after perform_all_checks() has run.
from naca_core.pre_run_checks import perform_all_checks

def plot_validation_results(results_dict, results_subfolder):
    """
    Saves two figures, one column per profile:
      validation_plot_cl.svg     Cl vs alpha: in-house, XFOIL inviscid (Mach 0), XFOIL viscous (real Mach)
      validation_plot_error.svg  error breakdown vs alpha: total = numerical + Mach + viscosity
    """
    import matplotlib.pyplot as plt

    os.makedirs(results_subfolder, exist_ok=True)
    n = len(results_dict)

    # --- Cl vs alpha ---
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), squeeze=False)
    for ax, (profile_name, data) in zip(axes[0], results_dict.items()):
        for key, style, label in (('Cl_InHouse', 'o--', 'In-house (potential)'),
                                  ('Cl_XFOIL_Inv_M0', 's-', 'XFOIL inviscid, Mach 0'),
                                  ('Cl_XFOIL_Visc_M', '^-', 'XFOIL viscous, real Mach')):
            pts = [(d['Alpha'], d[key]) for d in data if d[key] is not None]
            if pts:
                ax.plot(*zip(*pts), style, markersize=4, label=label)
        ax.set_title(profile_name)
        ax.set_xlabel('Alpha (degrees)')
        ax.set_ylabel('Cl')
        ax.grid(True)
        ax.legend(fontsize=8)
    fig.tight_layout()
    cl_file = os.path.join(results_subfolder, "validation_plot_cl.svg")
    fig.savefig(cl_file, format='svg')
    plt.close(fig)
    print(f"[i] Plot saved to: {cl_file}")

    # --- Error breakdown vs alpha ---
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), squeeze=False)
    for ax, (profile_name, data) in zip(axes[0], results_dict.items()):
        alphas = [d['Alpha'] for d in data]
        for key, marker, linestyle, color, label in (('dCl_Total', 'o', '-', 'black', 'Total'),
                                                     ('dCl_Numerical', 's', '--', 'tab:blue', 'Numerical'),
                                                     ('dCl_Mach', '^', '--', 'tab:orange', 'Mach'),
                                                     ('dCl_Viscosity', 'v', '--', 'tab:green', 'Viscosity')):
            # NaN for missing values: the line is interrupted instead of joining the points around a gap
            values = [d[key] if d[key] is not None else float('nan') for d in data]
            ax.plot(alphas, values, marker=marker, linestyle=linestyle, color=color, markersize=4, label=label)
            # Points computed with a single order (Mach first / visc first): hollow marker
            if key in ('dCl_Mach', 'dCl_Viscosity'):
                fallback = [d for d in data if d[key] is not None and d['Breakdown'] != 'average']
                if fallback:
                    ax.plot([d['Alpha'] for d in fallback], [d[key] for d in fallback], linestyle='none',
                            marker=marker, markersize=9, markerfacecolor='none', markeredgecolor=color)
        ax.axhline(0.0, color='grey', linewidth=0.8)
        ax.set_title(profile_name)
        ax.set_ylabel('dCl = in-house - reference')
        ax.set_xlabel('Alpha (degrees)  -  hollow marker: single-order breakdown')
        ax.grid(True)
        ax.legend(fontsize=8)
    fig.tight_layout()
    err_file = os.path.join(results_subfolder, "validation_plot_error.svg")
    fig.savefig(err_file, format='svg')
    plt.close(fig)
    print(f"[i] Plot saved to: {err_file}")

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
    import numpy as np

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
    perform_all_checks(require_xfoil=True)

    from naca_core.airfoil import naca4_airfoil
    from naca_core.panel_method import run_panel_analysis
    from naca_core.xfoil import XFoilAnalysis

    alphas, reynolds, mach, num_points = get_user_input()

    # TEST MATRIX
    # We test some common NACA profiles: symmetric, mildly cambered, and highly cambered
    profiles = [
        (0.0, 0.0, 0.12, "0012"),  # Symmetric
        (0.02, 0.4, 0.12, "2412"), # Standard cambered
        (0.04, 0.4, 0.12, "4412")  # Highly cambered
    ]

    # Four XFOIL runs per angle. Reference = XFOIL viscous at the real Mach number.
    xfoil_runs = {
        # key: (viscous, Mach, label used in the terminal)
        'inv_M0':  (False, 0.0,  "inviscid M0"),
        'inv_M':   (False, mach, "inviscid M"),
        'visc_M0': (True,  0.0,  "viscous M0"),
        'visc_M':  (True,  mach, "viscous M"),
    }

    # Lists and dictionaries to collect data
    all_results = []
    plot_data = {}

    n_cases = len(profiles) * len(alphas)
    print(f"Testing {len(profiles)} profiles and {len(alphas)} angles of attack "
          f"(Total cases: {n_cases}, 4 XFOIL runs each)")
    print("\n[i] dCl = in-house Cl minus XFOIL viscous Cl at the real Mach. Positive = the in-house solver overestimates Cl.")
    print("    total = num + Mach + visc")
    print("    num  = in-house vs XFOIL inviscid at Mach 0 (numerical error of the panel method)")
    print("    Mach = compressibility effect, visc = boundary-layer effect (average of the two orders)")
    print("    If an intermediate XFOIL run fails, a single order is used and the line is marked")
    print("    (Mach first / visc first / Mach+visc combined).")

    for m, p, t, name in profiles:
        print(f"\n--- Processing Profile NACA {name} ---")
        profile_name = f"NACA {name}"
        plot_data[profile_name] = []
        
        # 1. Geometry Generation
        X_panel, Y_panel, _ = naca4_airfoil(m, p, t, num_points=int(num_points/2)+1)
        
        for alpha in alphas:
            print(f"  > Alpha = {alpha:4.1f} deg ...", end=" ", flush=True)
            
            # --- In-House Solver (Panel Method) ---
            try:
                panel_res = run_panel_analysis(X_panel, Y_panel, float(alpha), verbose=False)
                cl_inhouse = panel_res['cl_potential']
            except Exception:
                cl_inhouse = None

            # --- XFOIL: four runs ---
            cl, cd_visc, problems = {}, None, []
            for key, (viscous, run_mach, label) in xfoil_runs.items():
                analysis = XFoilAnalysis(airfoil_name=f"naca_{name}", alpha=float(alpha),
                                         reynolds=reynolds, mach=run_mach, viscous=viscous)
                cl_run, cd_run, achieved_alpha = analysis.run_analysis(X_panel, Y_panel)
                if cl_run is None or abs(achieved_alpha - float(alpha)) > 0.1:
                    cl[key] = None
                    problems.append(f"{'TIMEOUT' if analysis.last_error == 'timeout' else 'FAILED'} {label}")
                else:
                    cl[key] = cl_run
                    if key == 'visc_M':
                        cd_visc = cd_run

            # --- Error breakdown ---
            # Best method allowed by the runs that succeeded; the parts always add up to the total.
            A, B, C, D, E = cl_inhouse, cl['inv_M0'], cl['inv_M'], cl['visc_M'], cl['visc_M0']
            d_total = A - D if None not in (A, D) else None
            d_num = d_mach = d_visc = d_mach_visc = None
            breakdown = "N/A"
            if None not in (A, B, D):
                d_num = A - B
                d_mach_visc = B - D  # Mach and viscosity together
                if C is not None and E is not None:
                    # Mach and viscosity interact: average of the two orders
                    d_mach = ((B - C) + (E - D)) / 2
                    d_visc = ((C - D) + (B - E)) / 2
                    breakdown = "average"
                elif C is not None:
                    d_mach, d_visc = B - C, C - D
                    breakdown = "Mach first"
                elif E is not None:
                    d_visc, d_mach = B - E, E - D
                    breakdown = "visc first"
                else:
                    breakdown = "Mach+visc combined"

            # Print quick result
            if cl_inhouse is None:
                problems.insert(0, "FAILED in-house")
            status = "OK" if not problems else ", ".join(problems)
            ih_str = f"{A:.4f}" if A is not None else "N/A"
            xf_str = f"{D:.4f}" if D is not None else "N/A"
            if breakdown in ("average", "Mach first", "visc first"):
                err_str = f"dCl total {d_total:+.4f} = num {d_num:+.4f} | Mach {d_mach:+.4f} | visc {d_visc:+.4f}"
                if breakdown != "average":
                    err_str += f" ({breakdown})"
            elif breakdown == "Mach+visc combined":
                err_str = f"dCl total {d_total:+.4f} = num {d_num:+.4f} | Mach+visc {d_mach_visc:+.4f} (Mach+visc combined)"
            elif d_total is not None:
                err_str = f"dCl total {d_total:+.4f} (breakdown N/A)"
            else:
                err_str = "dCl N/A"
            print(f"[{status}] In-House Cl: {ih_str} | XFOIL Cl: {xf_str} | {err_str}")
            
            res_dict = {
                'Profile': profile_name,
                'Alpha': float(alpha),
                'Cl_InHouse': A,
                'Cl_XFOIL_Inv_M0': B,
                'Cl_XFOIL_Inv_M': C,
                'Cl_XFOIL_Visc_M0': E,
                'Cl_XFOIL_Visc_M': D,
                'Cd_XFOIL_Visc_M': cd_visc,
                'dCl_Total': d_total,
                'dCl_Numerical': d_num,
                'dCl_Mach': d_mach,
                'dCl_Viscosity': d_visc,
                'dCl_Mach_plus_Viscosity': d_mach_visc,
                'Breakdown': breakdown,
            }
            
            all_results.append(res_dict)
            plot_data[profile_name].append(res_dict)

    # --- Save CSV ---
    # Same Results/ folder as the optimizer (repository root), in its own validation/ subfolder
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    results_subfolder_name = f"Re{int(round(reynolds))}_Mach{mach:.3f}"
    results_subfolder = os.path.join(repo_root, "Results", "validation", results_subfolder_name)
    os.makedirs(results_subfolder, exist_ok=True)
    csv_file = os.path.join(results_subfolder, "validation_results.csv")
    print(f"\n[+] Saving test results to '{csv_file}'...")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(all_results[0].keys()))
        writer.writeheader()
        writer.writerows(all_results)
        
    # --- Plot Generation ---
    print("[+] Generating comparative plots...")
    plot_validation_results(plot_data, results_subfolder)
        
    print("\n[+] Validation completed successfully!")

if __name__ == "__main__":
    main()
