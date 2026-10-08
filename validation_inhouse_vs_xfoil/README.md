# Validation: in-house vs XFOIL

Shows how far the in-house Cl is from XFOIL viscous, and **why**.

## Run

```bash
cd validation_inhouse_vs_xfoil
python run_validation.py
```

Press **Enter** to accept the default value of each prompt (fluid, angle-of-attack sweep, speed, chord, number of panels). XFOIL must be set up as described in `xfoil_viscous_optimizer/README.md`.

## What it computes

For NACA 0012, 2412 and 4412 at every angle, the script runs the in-house panel method (A) and XFOIL four times:

| | Run | Viscous | Mach |
|---|---|---|---|
| A | In-house | no | 0 |
| B | XFOIL | no | 0 |
| C | XFOIL | no | real |
| E | XFOIL | yes | 0 |
| D | XFOIL (reference) | yes | real |

The difference from the reference is split into three parts that always add up exactly to the total:

```
dCl total = A − D = num + Mach + visc
```

- **num** = A − B: numerical error of the panel method (same physics as XFOIL inviscid);
- **Mach**: compressibility, which the in-house solver ignores;
- **visc**: boundary layer, which the in-house solver ignores.

Mach and viscosity affect each other, so their parts are the average of the two possible orders. If an intermediate XFOIL run fails, the script uses the runs that worked and marks the line `(Mach first)`, `(visc first)` or `(Mach+visc combined)`; if XFOIL inviscid at Mach 0 fails it shows only the total, and if the reference run fails no error can be computed for that angle.

A positive dCl means the in-house solver overestimates Cl.

## Typical result

Default inputs (Re 3.4·10⁶, Mach 0.147, 160 panels), Windows, alpha = 8°:

| Airfoil | Total | Numerical | Mach | Viscosity |
|---|---|---|---|---|
| NACA 0012 | +0.054 | −0.008 | −0.015 | +0.077 |
| NACA 2412 | +0.083 | −0.012 | −0.018 | +0.112 |
| NACA 4412 | +0.127 | −0.016 | −0.020 | +0.162 |

The boundary layer explains most of the difference, and its share grows with camber and angle. The numerical error of the panel method is about 1% of Cl.

## Output

`Results/Validation_Re<Reynolds>_Mach<Mach>/`:

- `validation_results.csv`: all Cl values, the viscous Cd, the error parts and the breakdown method used;
- `validation_plot_cl.svg`: Cl vs alpha for each airfoil (in-house, XFOIL inviscid, XFOIL viscous);
- `validation_plot_error.svg`: the error parts vs alpha for each airfoil, with the total. Hollow markers = single-order breakdown.

## More

What every symbol, message and value means, and why each value was chosen: [Reference](../README.md#reference) in the main README.
