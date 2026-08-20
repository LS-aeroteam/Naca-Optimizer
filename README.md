# NACA Profile Optimizer

This project is a comprehensive suite for designing and optimizing 4-digit NACA airfoils. Born to meet various aerodynamic needs, the project is now organized into three main modules:

1. **inhouse_potential_optimizer**: The core of the project. It uses a potential flow solver (Source-Vortex Panel Method) developed entirely in-house to optimize and analyze airfoil geometry extremely fast, without any external dependencies.
2. **xfoil_viscous_optimizer**: The code that interfaces with the XFOIL executable to perform optimization taking into account viscosity and aerodynamic drag.
3. **validation_inhouse_vs_xfoil**: The environment dedicated to comparing and mathematically validating the results obtained by our in-house potential solver against XFOIL viscous data.
4. **_Original_projects**: Legacy code and original Matlab scripts (kept for reference).

Each folder contains its own README with usage instructions. To get started, we recommend trying out the inhouse_potential_optimizer!
