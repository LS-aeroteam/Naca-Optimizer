import matplotlib.pyplot as plt
import numpy as np

def plot_airfoil_geometry(X, Y, title="Airfoil Geometry", **kwargs):
    """
    Plots the geometry of the airfoil.

    Args:
        X (np.ndarray): X-coordinates of the airfoil.
        Y (np.ndarray): Y-coordinates of the airfoil.
        title (str): The title for the plot.
    """
    plt.figure(figsize=(12, 6))
    plt.plot(X, Y, 'k-', linewidth=1.5)
    plt.title(title, fontsize=16)
    plt.xlabel('x/c')
    plt.ylabel('y/c')
    plt.axis('equal')
    plt.grid(True)
    if 'save_path' in kwargs:
        plt.savefig(kwargs['save_path'], format='svg', bbox_inches='tight')
    plt.close()

def plot_pressure_coefficient(panel_results, alpha_deg, naca_name="", **kwargs):
    """
    Plots the pressure coefficient (Cp) distribution over the airfoil.

    Args:
        panel_results (dict): The results dictionary from `run_panel_analysis`.
        alpha_deg (float): The angle of attack for the title.
        naca_name (str): Optional name of the airfoil for the title.
    """
    Cp = panel_results['Cp']
    XC = panel_results['XC']
    num_panels = panel_results['num_panels']
    n_half = int(num_panels / 2)

    # Split into lower and upper surfaces (assuming clockwise point order)
    x_lower = XC[:n_half]
    cp_lower = Cp[:n_half]
    x_upper = XC[n_half:]
    cp_upper = Cp[n_half:]

    plt.figure(figsize=(10, 6))
    plt.plot(x_upper, cp_upper, 'b-o', markersize=3, linewidth=1.5, label='Upper Surface')
    plt.plot(x_lower, cp_lower, 'r-o', markersize=3, linewidth=1.5, label='Lower Surface')
    plt.gca().invert_yaxis()
    plt.xlabel('x/c')
    plt.ylabel('$C_p$')
    title = rf'Pressure Coefficient Distribution for {naca_name} $\alpha$ = {alpha_deg}° (Potential Analysis)'
    plt.title(title, fontsize=14)
    plt.grid(True)
    plt.legend()
    if 'save_path' in kwargs:
        plt.savefig(kwargs['save_path'], format='svg', bbox_inches='tight')
    plt.close()

def plot_lift_distribution(panel_results, alpha_deg, naca_name="", **kwargs):
    """
    Plots the lift distribution (Delta Cp) over the airfoil.

    Args:
        panel_results (dict): The results dictionary from `run_panel_analysis`.
        alpha_deg (float): The angle of attack for the title.
        naca_name (str): Optional name of the airfoil for the title.
    """
    Cp = panel_results['Cp']
    XC = panel_results['XC']
    num_panels = panel_results['num_panels']
    n_half = int(num_panels / 2)

    x_lower = XC[:n_half]
    cp_lower = Cp[:n_half]
    x_upper = XC[n_half:]
    cp_upper = Cp[n_half:]

    # x_upper goes from LE to TE. x_lower goes from TE to LE. 
    # Reverse lower to align with upper.
    x_aligned = x_upper
    delta_cp = cp_lower[::-1] - cp_upper

    plt.figure(figsize=(10, 6))
    plt.plot(x_aligned, delta_cp, 'g-o', markersize=3, linewidth=1.5, label=r'$\Delta C_p$ (Lift Distribution)')
    plt.xlabel('x/c')
    plt.ylabel(r'$\Delta C_p = C_{p,lower} - C_{p,upper}$')
    title = rf'Lift Distribution for {naca_name} $\alpha$ = {alpha_deg}° (Potential Analysis)'
    plt.title(title, fontsize=14)
    plt.grid(True)
    plt.legend()
    if 'save_path' in kwargs:
        plt.savefig(kwargs['save_path'], format='svg', bbox_inches='tight')
    plt.close()

def plot_optimization_history(history, **kwargs):
    """
    Plots the evolution of parameters and score during optimization.

    Args:
        history (list): A list of lists, where each inner list contains
                        [eval_count, m, p, t, cl, cd, score].
    """
    if not history:
        print("History is empty, cannot generate plot.")
        return

    # Filter out rows with non-numeric data (e.g., 'Failed', 'Bounds') and convert to numpy array
    history_np = np.array([row for row in history if isinstance(row[4], (int, float))], dtype=float)
    if history_np.shape[0] == 0:
        print("No successful evaluations in history, cannot generate plot.")
        return
        
    eval_count = history_np[:, 0]
    m = history_np[:, 1]
    p = history_np[:, 2]
    t = history_np[:, 3]
    scores = history_np[:, 6]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    # Plot 1: Parameters m, p, t
    ax1.plot(eval_count, m, 'r-o', markersize=3, label='m (Max Camber)')
    ax1.plot(eval_count, p, 'g-o', markersize=3, label='p (Camber Position)')
    ax1.plot(eval_count, t, 'b-o', markersize=3, label='t (Thickness)')
    ax1.set_ylabel('Parameter Value')
    ax1.set_title('Optimization History: Parameters & Score', fontsize=16)
    ax1.grid(True)
    ax1.legend()

    # Plot 2: Objective Function Score
    ax2.plot(eval_count, scores, 'k-o', markersize=3, label='Objective Score')
    ax2.set_xlabel('Evaluation Number')
    ax2.set_ylabel('Score (log scale)')
    ax2.set_yscale('log')
    ax2.grid(True)
    ax2.legend()
    
    plt.tight_layout()
    if 'save_path' in kwargs:
        plt.savefig(kwargs['save_path'], format='png', bbox_inches='tight')
    plt.close()
