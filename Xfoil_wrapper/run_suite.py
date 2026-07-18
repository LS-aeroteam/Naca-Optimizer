import os
import csv
import sys
from naca_aero_suite.pre_run_checks import perform_all_checks
from naca_aero_suite.optimizer import NacaOptimizer
from naca_aero_suite.airfoil import naca4_airfoil, save_airfoil_coordinates
from naca_aero_suite.panel_method import run_panel_analysis
from naca_aero_suite.plotting import (
    plot_airfoil_geometry, 
    plot_pressure_coefficient, 
    plot_lift_distribution,
    plot_optimization_history
)

def get_fluid_selection():
    fluids = {
        '1': {'name': 'Aria (Standard SL)', 'viscosity': 1.46e-5, 'speed_of_sound': 340.3},
        '2': {'name': 'Acqua (20°C)', 'viscosity': 1.00e-6, 'speed_of_sound': 1482.0}
    }
    while True:
        print("\nSeleziona il fluido operativo:")
        print("1. Aria (Viscosità cinematica: 1.46e-5 m^2/s, Velocità del suono: 340.3 m/s)")
        print("2. Acqua (Viscosità cinematica: 1.00e-6 m^2/s, Velocità del suono: 1482.0 m/s)")
        choice = input("Scelta (1 o 2) [default: 1]: ").strip() or '1'
        if choice in fluids:
            return fluids[choice]
        print("[!] Errore: selezione non valida.")

def get_user_input():
    """Gets all necessary inputs from the user."""
    print("\n======================================================================")
    print("                SETUP FOR DESIGN & ANALYSIS")
    print("======================================================================")
    
    fluid = get_fluid_selection()
    
    try:
        speed = float(input("\nEnter design speed (m/s) [e.g., 50]: ") or 50)
        
        # Bounding box constraints
        print("\n--- Bounding Box Constraints ---")
        max_length = float(input("Enter maximum physical length (L_max) in meters [e.g., 1.0]: ") or 1.0)
        max_height = float(input("Enter maximum physical height (H_max) in meters [e.g., 0.2]: ") or 0.2)
        
        chord = float(input(f"Enter airfoil chord (m) (must be <= {max_length}) [e.g., {max_length}]: ") or max_length)
        if chord > max_length:
            print(f"[!] Warning: chord {chord} exceeds max length {max_length}. Constraining chord to {max_length}.")
            chord = max_length
            
        target_reynolds = (speed * chord) / fluid['viscosity']
        target_mach = speed / fluid['speed_of_sound']
        
        print(f"\n [i] Calculated Conditions ({fluid['name']}):")
        print(f"    - Reynolds Number: {target_reynolds:,.0f}")
        print(f"    - Mach Number: {target_mach:.3f}")
        print(f"    - Max Allowed Profile Height: {max_height} m")

        target_alpha = float(input("\nEnter target angle of attack (degrees) [e.g., 4.0]: ") or 4.0)
        target_cl = float(input("Enter target lift coefficient (Cl) [e.g., 0.8]: ") or 0.8)
        max_cd = float(input("Enter maximum tolerable drag coefficient (Cd) [e.g., 0.02]: ") or 0.02)
        
        return target_reynolds, target_mach, target_alpha, target_cl, max_cd, chord, max_height
    except ValueError:
        print("\n[!] Invalid input. Please enter numerical values.")
        return None

def main():
    """Main execution block for the aerodynamic suite."""
    # Run all dependency and environment checks first
    perform_all_checks()
    
    # Get aerodynamic targets from the user
    inputs = get_user_input()
    if inputs is None:
        sys.exit(1)

    target_reynolds, target_mach, target_alpha, target_cl, max_cd, chord, max_height = inputs
    
    # --- FASE 1: RICERCA PROFILO ---
    print("\n======================================================================")
    print("                FASE 1: RICERCA PROFILO")
    print("======================================================================")
    initial_guess = [0.02, 0.4, 0.12]  # Start with a NACA 2412
    bounds = [(0.0, 0.09), (0.1, 0.7), (0.05, 0.25)] # Sensible bounds for NACA 4-digits

    optimizer = NacaOptimizer(
        reynolds=target_reynolds,
        alpha=target_alpha,
        target_cl=target_cl,
        max_cd=max_cd,
        mach=target_mach,
        chord=chord,
        max_height=max_height
    )

    result = optimizer.optimize(initial_guess, bounds, max_iter=50)

    if not result.success and result.nfev == 0:
        print("\n[!] Optimization failed to start. Please check your XFOIL setup and input parameters.")
        return
        
    # --- Process and Save Optimization Results ---
    opt_m, opt_p, opt_t = result.x
    naca_opt_str = f"{int(round(opt_m*100))}{int(round(opt_p*10))}{int(round(opt_t*100)):02d}"
    print(f"\n[+] OPTIMIZATION COMPLETE. Best profile found: NACA {naca_opt_str}")

    results_dir = f"Results_Re{int(target_reynolds)}_Alpha{target_alpha}_Cl{target_cl}"
    os.makedirs(results_dir, exist_ok=True)
    print(f"[i] Saving results to '{results_dir}/'")

    # Save optimized airfoil coordinates
    X_opt, Y_opt, _ = naca4_airfoil(opt_m, opt_p, opt_t, num_points=200)
    airfoil_filename = os.path.join(results_dir, f"airfoil_NACA_{naca_opt_str}.dat")
    save_airfoil_coordinates(X_opt, Y_opt, airfoil_filename)

    # Save optimization history
    history = optimizer.get_optimization_history()
    history_filename = os.path.join(results_dir, "optimization_history.csv")
    with open(history_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Eval', 'm', 'p', 't', 'Cl', 'Cd', 'Score'])
        writer.writerows(history)
    
    # Plot and save the optimization history
    history_plot_filename = os.path.join(results_dir, "optimization_history.svg")
    plot_optimization_history(history, save_path=history_plot_filename)
    
    # --- FASE 2: ANALISI PROFILO ---
    print("\n======================================================================")
    print("           FASE 2: ANALISI PROFILO")
    print("======================================================================")
    
    run_post = input("Run advanced analysis and generate plots? (y/n): ")
    if run_post.lower().strip() != 'y':
        return

    # 1. Plot airfoil geometry
    airfoil_plot_filename = os.path.join(results_dir, f"geometry_NACA_{naca_opt_str}.svg")
    plot_airfoil_geometry(X_opt, Y_opt, title=f"Optimized Airfoil: NACA {naca_opt_str}", save_path=airfoil_plot_filename)
    
    # --- FASE 3: ANALISI POTENZIALE E PRESSIONI ---
    print("\n======================================================================")
    print("            FASE 3: ANALISI POTENZIALE E PRESSIONI")
    print("======================================================================")
    try:
        num_panels = int(input("Enter number of panels for analysis [e.g., 160]: ") or 160)
        
        # Regenerate airfoil with panel-specific points if needed
        X_panel, Y_panel, _ = naca4_airfoil(opt_m, opt_p, opt_t, num_points=int(num_panels/2)+1)
        
        panel_results = run_panel_analysis(X_panel, Y_panel, target_alpha)
        
        # Add final Cl/Cd from optimizer to the results for plotting
        final_opt_perf = next((h for h in reversed(history) if isinstance(h[4], float)), None)
        final_cl = final_opt_perf[4] if final_opt_perf else None
        final_cd = final_opt_perf[5] if final_opt_perf else None

        print("\n--- Risultati Globali ---")
        if final_cl is not None and final_cd is not None:
            print(f"XFOIL Viscous Cl: {final_cl:.4f}")
            print(f"XFOIL Viscous Cd: {final_cd:.5f}")
        print(f"Potential Cl:     {panel_results['cl_potential']:.4f}")
        print("-------------------------\n")

        cp_plot_filename = os.path.join(results_dir, f"pressure_distribution_NACA_{naca_opt_str}.svg")
        plot_pressure_coefficient(panel_results, target_alpha, f"NACA {naca_opt_str}", save_path=cp_plot_filename)
        
        lift_dist_filename = os.path.join(results_dir, f"lift_distribution_NACA_{naca_opt_str}.svg")
        plot_lift_distribution(panel_results, target_alpha, f"NACA {naca_opt_str}", save_path=lift_dist_filename)
        
        # Save CSV with Cp, Cl, Cd
        csv_filename = os.path.join(results_dir, f"aerodynamic_data_NACA_{naca_opt_str}.csv")
        with open(csv_filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["# Global Coefficients"])
            writer.writerow(["XFOIL_Cl", "XFOIL_Cd", "Potential_Cl"])
            writer.writerow([final_cl, final_cd, panel_results['cl_potential']])
            writer.writerow([])
            writer.writerow(["# Surface Distributions"])
            writer.writerow(["x/c", "Cp_Upper", "Cp_Lower", "Delta_Cp"])
            
            n_half = panel_results['num_panels'] // 2
            x_upper = panel_results['XC'][n_half:]
            cp_upper = panel_results['Cp'][n_half:]
            cp_lower = panel_results['Cp'][:n_half][::-1]
            delta_cp = cp_lower - cp_upper
            
            for i in range(len(x_upper)):
                writer.writerow([f"{x_upper[i]:.6f}", f"{cp_upper[i]:.6f}", f"{cp_lower[i]:.6f}", f"{delta_cp[i]:.6f}"])

        print(f"\n[+] Analysis complete. All plots and CSV saved in '{results_dir}/'")

    except Exception as e:
        import traceback
        print(f"[!] An error occurred during panel analysis:")
        traceback.print_exc()

if __name__ == "__main__":
    main()

