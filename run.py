import os
import csv
import sys
import random

# Only the standard library is imported here: numpy/scipy/matplotlib are imported
# inside main(), after perform_all_checks() has verified (or installed) them.
from naca_core.pre_run_checks import perform_all_checks

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))

INHOUSE = "inhouse"
XFOIL = "xfoil"


def get_solver_selection():
    """Asks which solver to use. Returns INHOUSE or XFOIL."""
    solvers = {'1': INHOUSE, '2': XFOIL}
    while True:
        print("\nSelect the solver:")
        print("1. In-house panel method (potential flow, fast, no drag)")
        print("2. XFOIL (viscous, slower, gives drag)")
        choice = input("Choice (1 or 2) [default: 1]: ").strip()
        if not choice:
            return INHOUSE
        if choice in solvers:
            return solvers[choice]
        print("Error: invalid selection.")


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


def get_user_input(solver):
    """Gets all necessary inputs from the user. Returns a dict, or None on invalid input."""
    print("\n======================================================================")
    print("                SETUP FOR DESIGN & ANALYSIS")
    print("======================================================================")

    fluid = get_fluid_selection()
    fluid_viscosity = fluid['viscosity']
    speed_of_sound = fluid['speed_of_sound']

    try:
        speed = float(input("Enter design speed (m/s) [e.g., 50]: ") or 50)
        chord = float(input("Enter airfoil chord (m) [e.g., 1.0]: ") or 1.0)

        target_reynolds = (speed * chord) / fluid_viscosity
        target_mach = speed / speed_of_sound
        print(f"\n [i] Calculated Conditions:")
        print(f"    - Reynolds Number: {target_reynolds:,.0f}")
        print(f"    - Mach Number: {target_mach:.3f}")

        target_alpha = float(input("\nEnter target angle of attack (degrees) [e.g., 4.0]: ") or 4.0)

        # The objective can be chosen only with a solver that computes Cd (XFOIL for now)
        objective_choice, target_cl, max_cd = "1", None, None
        if solver == XFOIL:
            print("\nSelect the optimization objective:")
            print("1. Minimum Cd at a target Cl")
            print("2. Maximum Cl with a Cd limit")
            objective_choice = input("Choice (1 or 2) [default: 1]: ").strip() or "1"
            if objective_choice not in ("1", "2"):
                raise ValueError("invalid objective")
        if objective_choice == "1":
            target_cl = float(input("Enter target lift coefficient (Cl) [e.g., 0.8]: ") or 0.8)
        else:
            max_cd = float(input("Enter maximum drag coefficient (Cd) [e.g., 0.02]: ") or 0.02)

        max_height_box = float(input("Enter Bounding Box max height (m) [e.g., 0.3]: ") or 0.3)
        num_panels = int(input("Enter number of panels (Low=60, Medium=100, High=160) [default: 160]: ") or 160)
        seed_text = input("Enter random seed [default: random]: ").strip()
        seed = int(seed_text) if seed_text else random.SystemRandom().randrange(1, 1_000_000)
        span_text = input(f"Enter span for the 3D export (m) [default: {chord} = chord]: ").strip()
        span = float(span_text) if span_text else chord

        return {
            'reynolds': target_reynolds, 'mach': target_mach, 'alpha': target_alpha,
            'objective_choice': objective_choice, 'target_cl': target_cl, 'max_cd': max_cd,
            'max_height_box': max_height_box, 'chord': chord, 'num_panels': num_panels,
            'seed': seed, 'span': span,
        }

    except ValueError:
        print("\n[!] Invalid input. Please enter numerical values.")
        return None


def main():
    """Main execution block for the aerodynamic suite."""
    # The solver is chosen first: the XFOIL check below is needed only for XFOIL
    solver = get_solver_selection()

    # Run all dependency and environment checks first
    perform_all_checks(require_xfoil=(solver == XFOIL))

    from naca_core.airfoil import naca4_airfoil, save_airfoil_coordinates
    from naca_core.optimization import print_phase
    from naca_core.export import export_all
    from naca_core.panel_method import run_panel_analysis, surface_distributions
    from naca_core.plotting import (
        plot_airfoil_geometry,
        plot_pressure_coefficient,
        plot_lift_distribution,
        plot_optimization_history
    )
    if solver == XFOIL:
        from naca_core.xfoil_optimizer import NacaOptimizer, MIN_CD, MAX_CL, CL_TOLERANCE
        from naca_core.xfoil import XFoilAnalysis
    else:
        from naca_core.inhouse_optimizer import NacaOptimizer

    # Get aerodynamic targets from the user
    inputs = get_user_input(solver)
    if inputs is None:
        sys.exit(1)

    target_reynolds = inputs['reynolds']
    target_mach = inputs['mach']
    target_alpha = inputs['alpha']
    target_cl = inputs['target_cl']
    max_cd = inputs['max_cd']
    max_height_box = inputs['max_height_box']
    chord = inputs['chord']
    num_panels = inputs['num_panels']
    seed = inputs['seed']
    span = inputs['span']

    # --- PHASE 1: AIRFOIL OPTIMIZATION ---
    print_phase("PHASE 1: AIRFOIL OPTIMIZATION")

    if solver == XFOIL:
        objective = MIN_CD if inputs['objective_choice'] == "1" else MAX_CL
        if objective == MIN_CD:
            print(f"[i] Objective: minimum Cd with Cl = {target_cl} +/- {CL_TOLERANCE} (score = Cd in counts + penalty)")
        else:
            print(f"[i] Objective: maximum Cl with Cd <= {max_cd} (score = -Cl + penalty)")
        optimizer = NacaOptimizer(
            reynolds=target_reynolds,
            alpha=target_alpha,
            objective=objective,
            target_cl=target_cl,
            max_cd=max_cd,
            mach=target_mach,
            chord=chord,
            max_height_box=max_height_box,
            num_panels=num_panels,
            seed=seed
        )
    else:
        optimizer = NacaOptimizer(
            reynolds=target_reynolds,
            alpha=target_alpha,
            target_cl=target_cl,
            max_height_box=max_height_box,
            chord=chord,
            num_panels=num_panels,
            seed=seed
        )

    result = optimizer.optimize(max_iter=50)

    # Abort if no airfoil was analysed successfully: result.x would only be the minimum of the penalties
    history = optimizer.get_optimization_history()
    if not any(isinstance(row[4], float) for row in history):
        print("\n[!] No airfoil could be analysed successfully, so there is no valid result.")
        print("    Check the solver setup and the input values, then run again.")
        return

    # --- Process and Save Optimization Results ---
    opt_m, opt_p, opt_t = result.x
    naca_opt_str = f"{int(round(opt_m*100))}{int(round(opt_p*10))}{int(round(opt_t*100)):02d}"
    print(f"\n[+] OPTIMIZATION COMPLETE in {optimizer.eval_count} iterations. Best profile found: NACA {naca_opt_str}")
    print(f"[i] Exact parameters: m = {opt_m:.4f}, p = {opt_p:.4f}, t = {opt_t:.4f}")
    print(f"[i] Random seed: {seed} (enter it at the seed prompt to repeat this run)")

    # One Results/ folder for both solvers: Re<Re>_Alpha<alpha>_<objective>_<solver>_seed<N>
    objective_tag = f"Cl{target_cl}" if target_cl is not None else f"CdMax{max_cd}"
    results_dir_name = f"Re{int(round(target_reynolds))}_Alpha{target_alpha}_{objective_tag}_{solver}_seed{seed}"
    results_dir = os.path.join(REPO_ROOT, "Results", results_dir_name)
    os.makedirs(results_dir, exist_ok=True)
    print(f"[i] Saving results to 'Results/{results_dir_name}/'")

    # Save optimized airfoil coordinates
    X_opt, Y_opt, _ = naca4_airfoil(opt_m, opt_p, opt_t, num_points=200)
    airfoil_filename = os.path.join(results_dir, f"airfoil_NACA_{naca_opt_str}.dat")
    save_airfoil_coordinates(X_opt, Y_opt, airfoil_filename,
                             header=f"NACA {naca_opt_str} (m={opt_m:.6f} p={opt_p:.6f} t={opt_t:.6f})")

    # Save optimization history
    history_filename = os.path.join(results_dir, "optimization_history.csv")
    with open(history_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Eval', 'm', 'p', 't', 'Cl', 'Cd', 'Score'])
        writer.writerows(history)

    # Plot and save the optimization history
    history_plot_filename = os.path.join(results_dir, "optimization_history.svg")
    plot_optimization_history(history, save_path=history_plot_filename)

    # --- PHASE 2: AIRFOIL ANALYSIS ---
    print_phase("PHASE 2: AIRFOIL ANALYSIS")

    print("\n[+] Generating plots and extracting final aerodynamic data...")

    # 1. Plot airfoil geometry
    airfoil_plot_filename = os.path.join(results_dir, f"geometry_NACA_{naca_opt_str}.svg")
    plot_airfoil_geometry(X_opt, Y_opt, title=f"Optimized Airfoil: NACA {naca_opt_str}", save_path=airfoil_plot_filename)

    # --- PHASE 3: POTENTIAL & PRESSURE ANALYSIS ---
    print_phase("PHASE 3: POTENTIAL & PRESSURE ANALYSIS")
    try:
        # Regenerate airfoil with panel-specific points if needed
        X_panel, Y_panel, _ = naca4_airfoil(opt_m, opt_p, opt_t, num_points=int(num_panels/2)+1)

        panel_results = run_panel_analysis(X_panel, Y_panel, target_alpha)

        if solver == XFOIL:
            # Final XFOIL run on the optimized airfoil. The last row of the history is not used:
            # it can be a finite-difference or line-search point instead of result.x.
            # Same geometry (and paneling) as in the optimizer
            X_xfoil, Y_xfoil, _ = naca4_airfoil(opt_m, opt_p, opt_t, num_points=int(num_panels/2)+1)
            final_analysis = XFoilAnalysis(airfoil_name="final", alpha=target_alpha,
                                           reynolds=target_reynolds, mach=target_mach)
            final_cl, final_cd, final_alpha = final_analysis.run_analysis(X_xfoil, Y_xfoil)

        print("\n--- Global Results ---")
        if solver == XFOIL:
            if final_cl is not None and final_cd is not None:
                print(f"XFOIL Viscous Cl: {final_cl:.4f}")
                print(f"XFOIL Viscous Cd: {final_cd:.5f}")
                if abs(final_alpha - target_alpha) > 1e-6:
                    print(f"[!] XFOIL converged only up to alpha = {final_alpha} deg (target {target_alpha} deg)")
                if objective == MIN_CD and abs(final_cl - target_cl) > CL_TOLERANCE:
                    print(f"[!] Final Cl is outside the tolerance ({target_cl} +/- {CL_TOLERANCE})")
                if objective == MAX_CL and final_cd > max_cd:
                    print(f"[!] Final Cd is above the limit ({max_cd})")
            else:
                print("[!] XFOIL did not converge on the final airfoil: viscous Cl and Cd not available")
            print(f"Potential Cl:     {panel_results['cl_potential']:.4f}")
        else:
            print(f"Potential Cl:     {panel_results['cl_potential']:.4f}")
            print("Cd:               not available (potential flow)")
        print("-------------------------\n")

        cp_plot_filename = os.path.join(results_dir, f"pressure_distribution_NACA_{naca_opt_str}.svg")
        plot_pressure_coefficient(panel_results, target_alpha, f"NACA {naca_opt_str}", save_path=cp_plot_filename)

        lift_dist_filename = os.path.join(results_dir, f"lift_distribution_NACA_{naca_opt_str}.svg")
        plot_lift_distribution(panel_results, target_alpha, f"NACA {naca_opt_str}", save_path=lift_dist_filename)

        # Save CSV with Cp, Cl (and XFOIL Cl, Cd)
        csv_filename = os.path.join(results_dir, f"aerodynamic_data_NACA_{naca_opt_str}.csv")
        with open(csv_filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["# Global Coefficients"])
            if solver == XFOIL:
                writer.writerow(["XFOIL_Cl", "XFOIL_Cd", "Potential_Cl", "m", "p", "t", "Seed"])
                writer.writerow([final_cl, final_cd, panel_results['cl_potential'], opt_m, opt_p, opt_t, seed])
            else:
                writer.writerow(["Potential_Cl", "m", "p", "t", "Seed"])
                writer.writerow([panel_results['cl_potential'], opt_m, opt_p, opt_t, seed])
            writer.writerow([])
            writer.writerow(["# Surface Distributions (Cp_Lower interpolated at the upper-surface x/c)"])
            writer.writerow(["x/c", "Cp_Upper", "Cp_Lower", "Delta_Cp"])

            x_upper, cp_upper, cp_lower, delta_cp = surface_distributions(panel_results)

            for i in range(len(x_upper)):
                writer.writerow([f"{x_upper[i]:.6f}", f"{cp_upper[i]:.6f}", f"{cp_lower[i]:.6f}", f"{delta_cp[i]:.6f}"])

        # --- CAD / ParaView / OpenFOAM export ---
        print("\n[+] Exporting files for CAD, ParaView and OpenFOAM...")
        export_paths = export_all(X_opt, Y_opt, panel_results, chord, span, results_dir,
                                  name=f"airfoil_NACA_{naca_opt_str}")
        print(f"    - CAD (mm):        {os.path.basename(export_paths['dxf'])}, {os.path.basename(export_paths['csv'])}")
        print(f"    - ParaView (m):    {os.path.basename(export_paths['vtk'])}")
        print(f"    - OpenFOAM/3D (m): {os.path.basename(export_paths['stl'])} (span {span} m)")

        print(f"\n[+] Analysis complete. All plots and CSV saved in '{results_dir}/'")

    except Exception as e:
        import traceback
        print(f"[!] An error occurred during panel analysis:")
        traceback.print_exc()


if __name__ == "__main__":
    main()
