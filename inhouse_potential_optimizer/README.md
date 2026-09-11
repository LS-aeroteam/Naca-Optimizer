# In-house potential optimizer

Finds a NACA 4-digit airfoil that reaches a target lift coefficient, using our own source + vortex panel method (potential flow). No external programs are needed: only Python with `numpy`, `scipy` and `matplotlib`.

## Run

```bash
cd inhouse_potential_optimizer
python run.py
```

Press **Enter** to accept the default value of each prompt. If a Python library is missing, the script asks whether to install it with pip.

## What it does

1. Asks for fluid, speed, chord, angle of attack, target Cl, bounding-box height, number of panels and random seed.
2. Searches `m`, `p` and `t` in two phases:
   - **Genetic Algorithm** (scipy differential evolution): a random population of airfoils, always including a NACA 2412, explores the whole search range;
   - **SLSQP**: a gradient-based method refines the best airfoil found.

   Airfoils taller than the bounding box are not analysed (`OUT` in the table).
3. Analyses the best airfoil again and saves geometry, plots and data.

With the default inputs a full run takes a few seconds (about 3 s on our Linux test machine).

The same seed always gives the same run: the seed used is printed at the end and saved in the CSV.

## Limits

- **Potential flow:** no drag, no stall, no boundary layer. Reynolds and Mach are shown but not used. Cl is usually higher than the real (viscous) value: the validation script shows by how much and why.
- **Only Cl is matched:** many airfoils give the same Cl, so different seeds can return different airfoils. Minimising drag needs a boundary-layer model, which is planned.

## Output

`Results/Results_Re<Reynolds>_Alpha<angle>_Cl<target>/` with the airfoil coordinates (`.dat`), plots (`.svg`), the optimization history and the aerodynamic data (`.csv`). The full list of files is in the main README.

## More

What every symbol, message and value means, and why each value was chosen: [Reference](../README.md#reference) in the main README.
