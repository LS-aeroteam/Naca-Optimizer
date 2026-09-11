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
- **XFOIL:** viscous mode at your Reynolds and Mach number, with free transition (Ncrit = 9). XFOIL analyses the same points as the in-house solver (no re-paneling), so the number of panels you choose applies to both. To help convergence, the angle of attack is reached in 1° steps before the target value.

### 3. Optimization
The program changes `m`, `p` and `t` to make this score as small as possible:

```
score = (10 · (Cl − Cl_target))²  +  penalties
```

Penalties are added when the solver fails and (XFOIL only) when Cd is above your maximum. An airfoil taller than your bounding box is not analysed at all: it gets a penalty that grows with the excess height.

Both solvers use the same two-phase search:

1. **Genetic Algorithm** (global search, scipy differential evolution): a random population of airfoils, which always includes a NACA 2412, evolves for a few generations.
2. **SLSQP** (local refinement): a gradient-based method starts from the best airfoil of phase 1. The finite-difference step is 1e-4 for the in-house solver and 2e-3 for XFOIL, because XFOIL prints Cl with only 4 decimals.

The random population depends on a seed: the same seed gives exactly the same run.

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

Every script checks the Python libraries at start-up and, if some are missing, asks whether to install them with pip. The XFOIL solver and the validation script also check for XFOIL.

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
| Random seed | random | Enter a number to repeat a run exactly. The seed used is printed at the end and saved in the CSV |

While it runs, you see one table row per evaluated airfoil:

```
| Eval |   m    |   p    |   t    |   Cl    |   Cd    |  BB  |   Score    |
|    1 | 0.0200 | 0.4000 | 0.1200 |  0.7359 |       - |      | 4.1077e-01 |
|    2 | 0.0441 | 0.1438 | 0.1791 |  1.0147 |       - |      | 4.6090e+00 |
```

`BB` shows `OUT` when the airfoil is taller than the bounding box (it is not analysed). Cd shows `-` because potential flow has no drag.

### XFOIL solver

```bash
cd xfoil_viscous_optimizer
python run.py
```

Same prompts and same table as the in-house solver, plus **maximum drag coefficient** (default 0.02). In the Cl column, `Failed` means XFOIL did not converge at the target angle, and `Timeout` means it took more than 30 s. Both are treated as failed evaluations.

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
| `airfoil_NACA_xxxx.dat` | Airfoil coordinates (x y), XFOIL-compatible. The first line is the name with the exact m, p, t |
| `geometry_NACA_xxxx.svg` | Airfoil shape |
| `pressure_distribution_NACA_xxxx.svg` | Cp on upper and lower surface |
| `lift_distribution_NACA_xxxx.svg` | ΔCp = Cp,lower − Cp,upper along the chord (both taken at the same x/c) |
| `optimization_history.csv` | Every evaluated airfoil with its Cl, Cd and score |
| `optimization_history.svg` | How m, p, t and the score changed during the search |
| `aerodynamic_data_NACA_xxxx.csv` | Global coefficients, exact m, p, t and Cp distribution |

The validation script saves `validation_results.csv` and `validation_plot_cl.svg` in `Results/Validation_Re<Reynolds>_Mach<Mach>/`.

---

## Checking your changes

If you touch the geometry or the panel method, run:

```bash
python tests/baseline_check.py
```

It compares Cl and the full Cp distribution of NACA 0012, 2412 and 4412 at three angles of attack with a saved reference. If everything matches it ends with `All 24 quantities match the baseline`. If a change in results is intended, delete `tests/baseline_reference.json` and create a new one with `python tests/baseline_check.py --generate`.

---

## Limitations

- **Potential flow has no drag.** The in-house solver cannot predict drag, stall or the loss of lift caused by the boundary layer. Its Cl is usually higher than the viscous one.
- **Reynolds and Mach don't change the in-house result.** They are shown on screen and used in the folder name, but the panel method does not use them (no compressibility correction).
- **Many shapes give the same Cl.** The in-house search only matches the target Cl (plus the bounding box), so the airfoil it returns is one valid answer among many. The global search is random, so two runs can return different airfoils unless you enter the same seed.
- **The NACA name is rounded.** The optimizer works with continuous values (for example `m = 0.0190`, `p = 0.529`, `t = 0.1975` is saved as "NACA 2520"). The saved geometry uses the exact values, which are printed at the end of the run and written in the `.dat` header and in the CSV.
- **XFOIL does not always converge.** Airfoils where it fails get a large penalty and the search moves on.
- **Only NACA 4-digit airfoils.**

## Roadmap

- Boundary-layer model for the in-house solver, to estimate drag and transition without XFOIL
- One shared core package for both solvers

## Credits

- XFOIL is developed by Mark Drela and Harold Youngren (MIT) and released under the GNU General Public License.
