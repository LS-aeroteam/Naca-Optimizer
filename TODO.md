## Geometry (Geometry Engine)
- [ ]  **NACA 5-digit**: Implement the mathematical formulation that generates the NACA 5-digit series airfoils.
- [ ]  **Asymmetric airfoils**: Support independent curves for the upper and lower surface.
- [ ]  **Trailing-edge closure**: Set the geometric standard for the trailing edge (open vs closed), based on the requirements of the later mesh generation.

## Aerodynamics and Physical Models
- [ ]  **Boundary layer (Prandtl)**: Add the boundary-layer calculation to the in-house solver to estimate Cd. Goal: use the same cost function as XFOIL.
- [ ]  **Transition analysis (Spike)**: Study how XFOIL implements the transition criterion (e^N method or similar) and define the architecture to replicate/adapt it in the code.
- [ ]  **Separation handling (Spike)**: Design the algorithm that detects flow separation and accounts for its effects.
- [ ]  **Streamlines**: Output the streamlines for the winning airfoil.
- [ ]  **Prandtl-Glauert correction**: Add the correction for subsonic compressibility effects (baseline: the Mach number changes Cl by about -0.02 at 10°).

## I/O Architecture and Optimization
- [x]  **Objective 2 tuning**: Recalibrate the default Cd limit (the current threshold of 0.02 at 4° has no effect on the optimization).
- [x]  **Input file**: Define the inputs in a file instead of the terminal.
- [x]  **Solver selection from the CLI**: Evaluate adding a command-line option next to the interactive question (e.g. `--solver xfoil`), useful to run several cases automatically.
- [x]  **Terminal non-regression test**: Evaluate making the before/after comparison of the terminal output (fixed seed) a permanent test in `tests/`, next to the existing numerical baseline.
- [x]  **Cl penalty in objective 1**: With the quadratic penalty, going slightly outside the Cl tolerance costs less than the Cd saved (0.0003 outside the tolerance = 9 drag counts): with seed 657086 the final Cl is 0.7947 instead of 0.8 ± 0.005. Evaluate a linear penalty or a larger weight. It changes the results: to be decided together.
- [x]  **XFOIL logging**: `naca_core/xfoil.py` calls `logging.basicConfig(level=INFO)` when the module is loaded, so in XFOIL runs the `[INFO]` messages of ezdxf are shown during the export. Use a module logger without configuring the global one.

## Repository Management
- [ ]  **License**: Choose and apply an open-source license (e.g. MIT or GPL) by adding a `LICENSE` file, to define legally how third parties may use the code.
