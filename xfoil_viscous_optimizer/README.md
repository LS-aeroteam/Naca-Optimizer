# XFOIL viscous optimizer

Same workflow and same terminal output as the in-house optimizer, but every airfoil is analysed with [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) in viscous mode (your Reynolds and Mach number, free transition with Ncrit = 9). It gives both Cl and Cd.

## XFOIL setup

- **Windows:** the script looks for `xfoil.exe` in this folder. If it is not there, it downloads XFOIL 6.99 from the official MIT page.
- **Linux / macOS:** install XFOIL and make sure `xfoil` is on your PATH, or copy the executable into this folder.

## Run

```bash
cd xfoil_viscous_optimizer
python run.py
```

Press **Enter** to accept the default value of each prompt.

## Objectives

At start you choose one of two objectives:

1. **Minimum Cd at a target Cl.** Cl must stay within `target ± 0.005`; among those airfoils, the one with the lowest Cd wins.
2. **Maximum Cl with a Cd limit.** Among the airfoils with Cd below your limit, the one with the highest Cl wins.

The search is the same as the in-house optimizer: Genetic Algorithm, then SLSQP.

## How XFOIL is used

- XFOIL analyses the same points as the panel method (no re-paneling), so the number of panels you choose applies to both solvers.
- The target angle is reached in 1° steps. If XFOIL does not converge at the target angle, the evaluation counts as failed (`Failed` in the table). A run longer than 30 s is stopped (`Timeout`).
- XFOIL graphics are turned off: no plot windows open.
- The final Cl and Cd come from a new XFOIL run on the optimized airfoil.

## Speed

Every evaluation starts XFOIL, so a run is slower than the in-house one: with the default inputs, about 100 evaluations and 1–2 minutes (90 s measured on Windows).

## Output

`Results/Results_Re<Reynolds>_Alpha<angle>_Cl<target>/` (objective 1) or `.../Results_Re<Reynolds>_Alpha<angle>_CdMax<limit>/` (objective 2), with the same files as the in-house optimizer plus the XFOIL Cl and Cd in the aerodynamic CSV.

The `export/` subfolder contains the airfoil ready for CAD (DXF, CSV in mm), ParaView (VTK) and OpenFOAM (STL): see "Using the exported files" in the main README.

## More

What every symbol, message and value means, and why each value was chosen: [Reference](../README.md#reference) in the main README.
