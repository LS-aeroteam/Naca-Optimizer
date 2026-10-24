## Aerodynamics and Physical Models
- [ ]  **Boundary layer (Prandtl)**: Add the boundary-layer calculation to the in-house solver to estimate Cd. Goal: use the same cost function as XFOIL.
- [ ]  **Transition analysis (Spike)**: Study how XFOIL implements the transition criterion (e^N method or similar) and define the architecture to replicate/adapt it in the code.
- [ ]  **Separation handling (Spike)**: Design the algorithm that detects flow separation and accounts for its effects.
- [ ]  **Streamlines**: Output the streamlines for the winning airfoil.
- [ ]  **Prandtl-Glauert correction**: Add the correction for subsonic compressibility effects (baseline: the Mach number changes Cl by about -0.02 at 10°).

## Repository Management
- [ ]  **License**: Choose and apply an open-source license (e.g. MIT or GPL) by adding a `LICENSE` file, to define legally how third parties may use the code.
