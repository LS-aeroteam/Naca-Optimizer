# NACA Optimization Tool & Aerodynamic Suite

## 📖 The Project

This open-source repository contains a unified Python suite for the analysis, inverse design, and optimization of NACA 4-digit airfoils. Originally developed as separate projects, the code has been refactored into a single, modular, and robust tool that combines a viscous flow optimizer (using XFOIL) with an inviscid panel method for detailed aerodynamic analysis.

## ✨ Code Refactoring and Improvements

The original `NACA_aero_suite.py` script has been completely refactored into a structured Python package to improve maintainability, readability, and extensibility.

Key improvements include:
*   **Modular Structure:** The code is now organized into a `naca_aero_suite` package with clear separation of concerns (airfoil generation, XFOIL execution, optimization, panel method, and plotting).
*   **English Translation:** All code, comments, and user-facing prompts have been translated to English.
*   **Robust Error Handling:** The XFOIL wrapper now has improved error and timeout handling, preventing silent failures.
*   **Safe File Management:** The use of temporary directories for analysis files prevents race conditions and ensures that intermediate files are automatically cleaned up.
*   **Unified Geometry Generation:** A single, consistent function (`naca4_airfoil`) is now used for all airfoil geometry creation.
*   **New Features:** The suite now automatically saves optimization history to a `.csv` file and can plot the convergence of the optimizer.

## 📂 New Code Structure

The project is now organized as follows:

```
.
├── naca_aero_suite/
│   ├── __init__.py
│   ├── airfoil.py          # Generates NACA 4-digit airfoil coordinates.
│   ├── xfoil.py            # A robust wrapper for running XFOIL analysis.
│   ├── optimizer.py        # Manages the SLSQP optimization process.
│   ├── panel_method.py     # Implements the inviscid panel method.
│   └── plotting.py         # Contains all plotting functions.
├── run_suite.py            # Main user-facing script to run the application.
├── NACA_aero_suite.py      # The original, monolithic script (archived).
└── README.md
```

## 🛠️ Prerequisites & Installation

Ensure you have **Python 3.8+** installed on your system.

1.  **Install Python Dependencies:**
    ```bash
    pip install numpy scipy matplotlib
    ```

2.  **XFOIL Installation:**
    The optimizer depends on the XFOIL executable.
    *   **Windows:** The suite will automatically download and install `xfoil.exe` for you during the first run! No manual action is required.
    *   **Linux/WSL:** Ensure `xfoil` is installed and accessible from your system's PATH (e.g., `sudo apt-get install xfoil`), or manually place the executable in the project's root directory.

## 🚀 How to Run

All functionality is now accessed through the `run_suite.py` script.

1.  Open your terminal in the project's root directory.
2.  Run the script:
    ```bash
    python run_suite.py
    ```

The script will guide you through a command-line interface:

**Phase 1: Optimization**
You will be prompted to enter the design conditions and aerodynamic targets:
*   Design speed and chord (to calculate Reynolds and Mach numbers).
*   Target angle of attack ($\alpha$).
*   Target lift coefficient ($C_L$).
*   Maximum allowable drag coefficient ($C_D$).

The optimizer will then run to find the NACA 4-digit airfoil that best meets these targets.

**Phase 2: Analysis and Visualization**
Once the optimization is complete, the results (airfoil `.dat` file, optimization history `.csv`) will be saved to a new `Results_*` directory.

You will then be prompted to run an advanced analysis on the optimized profile, which includes:
*   A plot of the optimization history.
*   A plot of the final airfoil geometry.
*   An inviscid panel method analysis, with a plot of the pressure coefficient ($C_p$) distribution.
