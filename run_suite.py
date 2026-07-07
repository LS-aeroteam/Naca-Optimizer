import os
import csv
from naca_aero_suite.optimizer import NacaOptimizer
from naca_aero_suite.airfoil import naca4_airfoil, save_airfoil_coordinates
from naca_aero_suite.panel_method import run_panel_analysis
from naca_aero_suite.plotting import (
    plot_airfoil_geometry, 
    plot_pressure_coefficient, 
    plot_optimization_history
)

def get_user_input():
    """Gets all necessary inputs from the user."""
    print("======================================================================")
    print("      AERODYNAMIC SUITE: INVERSE DESIGN & INVISCID ANALYSIS")
    print("======================================================================")
    
    # For simplicity, we'll use air at sea level properties.
    # This could be expanded to a fluid selection menu.
    fluid_viscosity = 1.46e-5  # Kinematic viscosity of air at SL in m^2/s
    speed_of_sound = 340.3     # Speed of sound in m/s
    
    try:
        speed = float(input("Enter design speed (m/s) [e.g., 50]: ") or 50)
        chord = float(input("Enter airfoil chord (m) [e.g., 1.0]: ") or 1.0)
        
        target_reynolds = (speed * chord) / fluid_viscosity
        target_mach = speed / speed_of_sound
        print(f"
[i] Calculated Conditions:")
        print(f"    - Reynolds Number: {target_reynolds:,.0f}")
        print(f"    - Mach Number: {target_mach:.3f}")

        target_alpha = float(input("
Enter target angle of attack (degrees) [e.g., 4.0]: ") or 4.0)
        target_cl = float(input("Enter target lift coefficient (Cl) [e.g., 0.8]: ") or 0.8)
        max_cd = float(input("Enter maximum tolerable drag coefficient (Cd) [e.g., 0.02]: ") or 0.02)
        
        return target_reynolds, target_mach, target_alpha, target_cl, max_cd
    except ValueError:
        print("
[!] Invalid input. Please enter numerical values.")
        return None

def main():
    """Main execution block for the aerodynamic suite."""
    
    inputs = get_user_input()
    if inputs is None:
        return

    target_reynolds, target_mach, target_alpha, target_cl, max_cd = inputs
    
    # --- Run Optimization ---
    initial_guess = [0.02, 0.4, 0.12]  # Starting with a NACA 2412
    bounds = [(0.0, 0.09), (0.1, 0.7), (0.05, 0.25)] # Sensible bounds for NACA 4-digits

    optimizer = NacaOptimizer(
        reynolds=target_reynolds,
        alpha=target_alpha,
        target_cl=target_cl,
        max_cd=max_cd,
        mach=target_mach
    )

    result = optimizer.optimize(initial_guess, bounds, max_iter=50)

    if not result.success and result.nfev == 0:
        print("
[!] Optimization failed to start. Please check your setup and XFOIL installation.")
        return
        
    # --- Process and Save Results ---
    opt_m, opt_p, opt_t = result.x
    naca_opt_str = f"{int(round(opt_m*100))}{int(round(opt_p*10))}{int(round(opt_t*100)):02d}"
    print(f"
[+] OPTIMIZATION COMPLETE. Best profile found: NACA {naca_opt_str}")

    # Create a directory for the results
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
    
    # --- Post-Analysis and Plotting ---
    run_post = input("
Do you want to run advanced analysis and view plots for the optimized profile? (y/n): ")
    if run_post.lower().strip() != 'y':
        return

    # 1. Plot optimization history
    plot_optimization_history(history)
    
    # 2. Plot airfoil geometry
    plot_airfoil_geometry(X_opt, Y_opt, title=f"Optimized Airfoil: NACA {naca_opt_str}")
    
    # 3. Run Panel Method and Plot Cp
    try:
        num_panels = int(input("
Enter number of panels for analysis [e.g., 160]: ") or 160)
        
        # Panel method needs coordinates in clockwise order
        # Our naca4_airfoil function already provides this
        panel_results = run_panel_analysis(X_opt, Y_opt, target_alpha)
        
        plot_pressure_coefficient(panel_results, target_alpha, f"NACA {naca_opt_str}")

    except ValueError:
        print("[!] Invalid number of panels. Aborting analysis.")
    except Exception as e:
        print(f"[!] An error occurred during panel analysis: {e}")

if __name__ == "__main__":
    main()
