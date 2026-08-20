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
You will be prompted for a few straightforward parameters (fluid, speed, chord, target Cl, and maximum thickness). The software will run automatically, display the results in the terminal, and save the coordinates, vector plots, and textual data (CSV) in a dedicated folder.
