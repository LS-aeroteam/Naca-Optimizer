# NACA Profile Optimizer

This project is a comprehensive suite for designing and optimizing 4-digit NACA airfoils. Born to meet various aerodynamic needs, the project is now organized into three main modules:

1. **InHouse_Solver**: The core of the project. It uses a potential flow solver (Source-Vortex Panel Method) developed entirely in-house to optimize and analyze airfoil geometry extremely fast, without any external dependencies.
2. **Xfoil_wrapper**: The legacy code that interfaces with the XFOIL executable to perform optimization taking into account viscosity and aerodynamic drag.
3. **InHouse_Validation**: The environment dedicated to comparing and mathematically validating the results obtained by our In-House solver against XFOIL data or other references.

Each folder contains its own README with usage instructions. To get started, we recommend trying out the In-House Solver!
