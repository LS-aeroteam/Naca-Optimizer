# In-House Potential Optimizer

This module represents the mathematical core of the project. The script performs aerodynamic optimization using a **Source-Vortex Panel Method** originally developed in MATLAB and ported here to Python to achieve extremely high performance.

## Main Features
- **No External Dependencies**: It solves the potential flow field by calculating pressure (Cp) and lift (Cl) distributions completely autonomously, without needing to interface with XFOIL or other executables.
- **Hybrid Optimization**: The search for the ideal airfoil is a two-step process. It starts with a Genetic Algorithm to widely explore the design space and avoid getting stuck in local minima, followed by a sub-millimeter gradient-based refinement (SLSQP method).
- **Bounding Box**: You can impose a maximum wing height. The solver will penalize and discard geometries that, once scaled by the chord, turn out to be too thick for your project constraints.

## Usage
Simply run the main script from the terminal:
```bash
python run.py
```
Upon launching, the pre-run checks will verify your environment and **automatically download and install** any missing Python dependencies (like `numpy`, `scipy`, or `matplotlib`) in the background.

You will be prompted for a few straightforward parameters (fluid, speed, chord, target angle of attack, target Cl, and maximum thickness constraint). The software will run automatically without interruptions, displaying the optimization progress directly in the terminal.

### Outputs
To keep your workspace clean, all outputs are isolated. The generated airfoil coordinates, pressure vector plots, and detailed textual data (CSV) are automatically saved into a categorized subfolder within a dedicated `Results/` directory.
