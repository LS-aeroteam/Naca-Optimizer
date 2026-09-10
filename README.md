# NACA Optimizer

Find a NACA 4-digit airfoil that gives you the lift you need.

You tell the program the flow conditions (fluid, speed, chord), the angle of attack and the lift coefficient you want. It searches the NACA 4-digit family for a shape that hits that target, then saves the coordinates, the plots and the raw data.

The project comes with two solvers that follow the same workflow:

- **In-house potential solver** – our own panel method, written in pure Python. No external programs needed. It gives lift and pressure distribution, but it ignores viscosity, so it cannot compute drag.
- **XFOIL viscous solver** – runs [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) in the background. Slower, but it includes the boundary layer, so it also gives drag and can respect a maximum drag coefficient.

A third script compares the two, so you can see how far potential flow is from the viscous result.

---

## Project structure

```
Naca-optimizer/
├── naca_core/                      # Shared code: geometry, panel method, XFOIL wrapper, plots, checks
├── inhouse_potential_optimizer/    # Optimizer based on our panel method
│   ├── run.py                      # Start here
│   ├── inhouse_optimizer.py
│   └── Results/                    # Output folders (one per run)
├── xfoil_viscous_optimizer/        # Optimizer based on XFOIL
│   ├── run.py
│   ├── xfoil_optimizer.py
│   └── Results/
├── validation_inhouse_vs_xfoil/    # In-house vs XFOIL comparison
│   ├── run_validation.py
│   └── Results/
├── tests/                          # Numerical baseline check
└── _Original_projects/             # Original MATLAB script and first Python version (reference only)
```

---

## How it works

### 1. Geometry
Each airfoil is built from three numbers:

| Parameter | Meaning | Example (NACA 2412) |
|-----------|---------|---------------------|
| `m` | maximum camber (fraction of chord) | 0.02 |
| `p` | position of maximum camber (fraction of chord) | 0.4 |
| `t` | maximum thickness (fraction of chord) | 0.12 |

Points are placed with cosine spacing, so there are more of them near the leading and trailing edge, where the flow changes fastest.

### 2. Aerodynamics
- **In-house:** a source + vortex panel method. The surface is split into straight panels; each panel carries its own source strength and all panels share one vortex strength. The Kutta condition at the trailing edge closes the system. Lift comes from the total circulation (Kutta–Joukowski). This solver was first written in MATLAB (`_Original_projects/Original_Matlab_Project/Prova_finale.m`) and then ported to Python.
- **XFOIL:** viscous mode at your Reynolds and Mach number, with free transition (Ncrit = 9). To help convergence, the angle of attack is reached in 1° steps before the target value.

### 3. Optimization
The program changes `m`, `p` and `t` to make this score as small as possible:

```
score = (10 · (Cl − Cl_target))²  +  penalties
```

Penalties are added when the airfoil is taller than your bounding box, when the solver fails, and (XFOIL only) when Cd is above your maximum.

- **In-house:** a short differential evolution run (a genetic-style global search) followed by SLSQP (a gradient-based local refinement).
- **XFOIL:** SLSQP starting from a NACA 2412.

Search limits: `0 ≤ m ≤ 0.09`, `0.1 ≤ p ≤ 0.7`, `0.05 ≤ t ≤ 0.25`.

---

## Requirements

- Python 3.12 or newer (tested with 3.12 and 3.14)
- `numpy`, `scipy`, `matplotlib` (minimum versions in `requirements.txt`)
- XFOIL – only for the XFOIL solver and the validation script

## Installation

```bash
git clone https://github.com/LS-aeroteam/Naca-Optimizer.git
cd Naca-Optimizer
python -m pip install -r requirements.txt
```

**XFOIL setup**

- **Windows:** the script looks for `xfoil.exe` inside `xfoil_viscous_optimizer/`. If it is not there, it downloads XFOIL 6.99 from the official MIT page and puts it there.
- **Linux / macOS:** install XFOIL yourself and make sure the `xfoil` command is on your PATH, or copy the executable into `xfoil_viscous_optimizer/`.

Every script checks the Python libraries at start-up and tries to install any that are missing. The XFOIL solver and the validation script also check for XFOIL.

---

## Usage

Every prompt has a default value: just press **Enter** to accept it.

### In-house solver

```bash
cd inhouse_potential_optimizer
python run.py
```

| Prompt | Default | Notes |
|--------|---------|-------|
| Fluid | 1 (air) | 1 = air at sea level, 2 = water at 20 °C |
| Design speed | 50 m/s | |
| Chord | 1.0 m | Reynolds and Mach are computed from speed, chord and fluid |
| Angle of attack | 4.0° | |
| Target Cl | 0.8 | |
| Bounding box max height | 0.3 m | Maximum total height of the airfoil (thickness + camber) |
| Number of panels | 160 | 60 = fast, 160 = accurate |

While it runs, you see one table row per evaluated airfoil:

```
| Eval |   m    |   p    |   t    |   Cl   | BB OOB |   Score    |
|    1 | 0.0146 | 0.6218 | 0.0973 | 0.7047 |        | 9.0794e-01 |
|    2 | 0.0424 | 0.3603 | 0.1974 | 1.0729 |        | 7.4451e+00 |
```

`BB OOB` shows `H` when the airfoil is taller than the bounding box.

### XFOIL solver

```bash
cd xfoil_viscous_optimizer
python run.py
```

Same prompts as the in-house solver, plus **maximum drag coefficient** (default 0.02). The table shows Cd next to Cl. A `*` next to a value means XFOIL did not converge exactly at the target angle and the closest converged angle was used.

### Validation

```bash
cd validation_inhouse_vs_xfoil
python run_validation.py
```

It runs NACA 0012, 2412 and 4412 over an angle-of-attack sweep (default −4° to 10°, step 2°) with both solvers and prints the lift coefficient from each one.

Keep in mind what is being compared: the in-house solver is potential flow, XFOIL here is viscous. A gap of a few percent in Cl is expected, because the boundary layer makes the real airfoil behave as if it had less camber. The gap tells you how much viscosity matters at your conditions; it is not only a numerical error.

---

## Output

Each optimization run creates its own folder:

```
Results/Results_Re<Reynolds>_Alpha<angle>_Cl<target>/
```

| File | Content |
|------|---------|
| `airfoil_NACA_xxxx.dat` | Airfoil coordinates (x y) |
| `geometry_NACA_xxxx.svg` | Airfoil shape |
| `pressure_distribution_NACA_xxxx.svg` | Cp on upper and lower surface |
| `lift_distribution_NACA_xxxx.svg` | ΔCp = Cp,lower − Cp,upper along the chord |
| `optimization_history.csv` | Every evaluated airfoil with its Cl, Cd and score |
| `optimization_history.svg` | How m, p, t and the score changed during the search |
| `aerodynamic_data_NACA_xxxx.csv` | Global coefficients and Cp distribution |

The validation script saves `validation_results.csv` and `validation_plot_cl.svg` in `Results/Validation_Re<Reynolds>_Mach<Mach>/`.

---

## Checking your changes

If you touch the geometry or the panel method, run:

```bash
python tests/baseline_check.py
```

It compares Cl and the full Cp distribution of NACA 0012, 2412 and 4412 at three angles of attack with a saved reference. If everything matches it ends with `All 24 quantities match the baseline`. If a change in results is intended, delete `tests/baseline_reference.json` and create a new one with `python tests/baseline_check.py --generate`.

---

## Known issues

These problems are already identified and fixes are planned:

- `airfoil_NACA_xxxx.dat` is written on a single line (the line breaks are missing), so XFOIL and other tools can't read it yet.
- `optimization_history.svg` is actually a PNG image saved with the wrong extension. Rename it to `.png` to open it.

## Limitations

- **Potential flow has no drag.** The in-house solver cannot predict drag, stall or the loss of lift caused by the boundary layer. Its Cl is usually higher than the viscous one.
- **Reynolds and Mach don't change the in-house result.** They are shown on screen and used in the folder name, but the panel method does not use them (no compressibility correction).
- **Many shapes give the same Cl.** The in-house search only matches the target Cl (plus the bounding box), so the airfoil it returns is one valid answer among many. The global search is random, so two runs can return different airfoils.
- **The NACA name is rounded.** The optimizer works with continuous values (for example `m = 0.0190`, `p = 0.529`, `t = 0.1975` is saved as "NACA 2520"), while the saved geometry uses the exact values.
- **XFOIL does not always converge.** Airfoils where it fails get a large penalty and the search moves on.
- **Only NACA 4-digit airfoils.**

## Roadmap

- Boundary-layer model for the in-house solver, to estimate drag and transition without XFOIL
- One shared core package for both solvers

## Credits

- XFOIL is developed by Mark Drela and Harold Youngren (MIT) and released under the GNU General Public License.
