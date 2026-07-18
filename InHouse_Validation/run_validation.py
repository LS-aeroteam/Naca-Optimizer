import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

# Aggiungiamo il percorso del Xfoil_wrapper per poter importare i suoi moduli completi
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Xfoil_wrapper")))

from naca_aero_suite.airfoil import naca4_airfoil
from naca_aero_suite.panel_method import run_panel_analysis
from naca_aero_suite.xfoil import XFoilAnalysis
from naca_aero_suite.pre_run_checks import perform_all_checks

def plot_validation_results(results_dict, alphas):
    """Genera grafici comparativi Cl vs Alpha."""
    plt.figure(figsize=(10, 6))
    
    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    
    for idx, (profile_name, data) in enumerate(results_dict.items()):
        c = colors[idx % len(colors)]
        
        # Filtra i valori validi
        valid_alphas = [d['Alpha'] for d in data if d['Cl_Xfoil'] is not None]
        cl_inhouse = [d['Cl_InHouse'] for d in data if d['Cl_Xfoil'] is not None]
        cl_xfoil = [d['Cl_Xfoil'] for d in data if d['Cl_Xfoil'] is not None]
        
        if not valid_alphas:
            continue
            
        plt.plot(valid_alphas, cl_inhouse, marker='o', linestyle='--', color=c, label=f'{profile_name} (In-House)')
        plt.plot(valid_alphas, cl_xfoil, marker='s', linestyle='-', color=c, label=f'{profile_name} (XFOIL)')

    plt.title('Confronto Coefficiente di Portanza (Cl) vs Angolo di Attacco (Alpha)')
    plt.xlabel('Alpha (gradi)')
    plt.ylabel('Cl')
    plt.grid(True)
    plt.legend()
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    plot_filename = os.path.join(base_dir, "validation_plot_cl.svg")
    plt.savefig(plot_filename)
    print(f"[i] Grafico salvato in: {plot_filename}")

def main():
    print("=========================================================")
    print("       IN-HOUSE SOLVER vs XFOIL VALIDATION SCRIPT")
    print("=========================================================")
    
    # Esegue i controlli preliminari (assicura che XFOIL sia presente)
    perform_all_checks()

    # MATRICE DI TEST
    # Testiamo alcuni profili NACA comuni: simmetrico, mediamente asimmetrico e molto curvo
    profiles = [
        (0.0, 0.0, 0.12, "0012"),  # Simmetrico
        (0.02, 0.4, 0.12, "2412"), # Asimmetrico standard
        (0.04, 0.4, 0.12, "4412")  # Asimmetrico molto curvo
    ]
    
    # Range di angoli di attacco
    alphas = np.arange(-4, 12, 2)  # -4, -2, 0, 2, 4, 6, 8, 10
    
    # Condizioni di flusso per XFOIL
    reynolds = 1e6
    mach = 0.0
    num_points = 160 # Numero di punti/pannelli
    
    # Liste e dizionari per raccogliere i dati
    all_results = []
    plot_data = {}

    print(f"Test su {len(profiles)} profili e {len(alphas)} angoli di attacco (Totale run: {len(profiles)*len(alphas)})")
    
    for m, p, t, name in profiles:
        print(f"\n--- Elaborazione Profilo NACA {name} ---")
        profile_name = f"NACA {name}"
        plot_data[profile_name] = []
        
        # 1. Generazione Geometria
        X_panel, Y_panel, _ = naca4_airfoil(m, p, t, num_points=int(num_points/2)+1)
        
        # 2. Inizializzazione Analizzatore XFOIL
        xfoil_analyzer = XFoilAnalysis(airfoil_name=f"naca_{name}", alpha=0.0, reynolds=reynolds, mach=mach)
        
        for alpha in alphas:
            print(f"  > Alpha = {alpha:2d}° ...", end=" ", flush=True)
            
            # --- In-House Solver (Metodo dei Pannelli) ---
            # run_panel_analysis nel wrapper Xfoil non supporta verbose=False, stamperà a schermo
            try:
                panel_res = run_panel_analysis(X_panel, Y_panel, float(alpha))
                cl_inhouse = panel_res['cl_potential']
            except Exception as e:
                cl_inhouse = 0.0
                print(f"[ERR Pannelli: {e}]", end=" ")
            
            # --- XFOIL ---
            xfoil_analyzer.alpha = float(alpha)
            cl_xfoil, cd_xfoil, achieved_alpha = xfoil_analyzer.run_analysis(X_panel, Y_panel)
            
            # --- Calcolo Errore ---
            if cl_xfoil is not None and cl_xfoil != 0:
                err_cl = abs((cl_inhouse - cl_xfoil) / cl_xfoil) * 100
            else:
                err_cl = None
                
            # Stampa risultato rapido
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

    # --- Salvataggio CSV ---
    base_dir = os.path.dirname(os.path.abspath(__file__))
    csv_file = os.path.join(base_dir, "validation_results.csv")
    print(f"\n[+] Salvataggio risultati test in '{csv_file}'...")
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['Profile', 'Alpha', 'Cl_InHouse', 'Cl_Xfoil', 'Error_Cl_%'])
        writer.writeheader()
        writer.writerows(all_results)
        
    # --- Generazione Grafici ---
    print(f"[+] Generazione grafici comparativi...")
    plot_validation_results(plot_data, alphas)
        
    print("\n[+] Validazione Completata con Successo!")

if __name__ == "__main__":
    main()
