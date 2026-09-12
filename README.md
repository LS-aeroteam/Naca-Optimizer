# NACA Optimizer

Find a NACA 4-digit airfoil that gives you the lift you need.

You tell the program the flow conditions (fluid, speed, chord), the angle of attack and the lift coefficient you want. It searches the NACA 4-digit family for a shape that hits that target, then saves the coordinates, the plots and the raw data.

The project comes with two solvers that follow the same workflow:

- **In-house potential solver** – our own panel method, written in pure Python. No external programs needed. It gives lift and pressure distribution, but it ignores viscosity, so it cannot compute drag.
- **XFOIL viscous solver** – runs [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) in the background. Slower, but it includes the boundary layer, so it also gives drag: it can find the lowest drag at a target lift, or the highest lift under a drag limit.

A third script compares the two, so you can see how far potential flow is from the viscous result, and why.

What every symbol, message and number means, and why each value was chosen, is explained in [Reference](#reference).

---

## Project structure

```
Naca-optimizer/
├── naca_core/                      # Shared code: geometry, panel method, XFOIL wrapper, plots, checks
├── inhouse_potential_optimizer/    # Optimizer based on our panel method
│   ├── run.py                      # Start here
│   ├── inhouse_optimizer.py
│   └── Results/                    # Output folders, one per run (created by the script, not versioned)
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
The program changes `m`, `p` and `t` to make a score as small as possible. The score depends on the solver:

- **In-house** (no drag available): reach the target Cl.
  ```
  score = (10 · (Cl − Cl_target))²
  ```
- **XFOIL**, two objectives chosen at start:
  1. **Minimum Cd at a target Cl.** The score is Cd in drag counts (1 count = 0.0001). If Cl is outside `Cl_target ± 0.005`, a penalty is added that quickly outweighs any Cd difference.
  2. **Maximum Cl with a Cd limit.** The score is −Cl. If Cd is above the limit, a penalty is added.

A large penalty is added when the solver fails. An airfoil taller than your bounding box is not analysed at all: it gets a penalty that grows with the excess height.

Both solvers use the same two-phase search:

1. **Genetic Algorithm** (global search, scipy differential evolution): a random population of airfoils, which always includes a NACA 2412, evolves for a few generations.
2. **SLSQP** (local refinement): a gradient-based method starts from the best airfoil of phase 1. The finite-difference step is 1e-4 for the in-house solver and 2e-3 for XFOIL, because XFOIL prints Cl with only 4 decimals.

The random population depends on a seed: the same seed gives exactly the same run.

Search limits: `0 ≤ m ≤ 0.09`, `0.1 ≤ p ≤ 0.7`, `0.05 ≤ t ≤ 0.25`.

---

## Requirements

- Python 3.12 or newer (tested with 3.12 and 3.14)
- `numpy`, `scipy`, `matplotlib`, `ezdxf` (for the DXF export); minimum versions in `requirements.txt`
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
| Span for the 3D export | same as the chord | Length of the extruded wing in the STL file, in metres |

While it runs, you see one table row per evaluated airfoil:

```
| Eval |   m    |   p    |   t    |   Cl    |   Cd    |  BB  |    Score    |
|    1 | 0.0200 | 0.4000 | 0.1200 |  0.7359 |       - |      |  4.1077e-01 |
|    2 | 0.0441 | 0.1438 | 0.1791 |  1.0147 |       - |      |  4.6090e+00 |
```

`BB` shows `OUT` when the airfoil is taller than the bounding box (it is not analysed). Cd shows `-` because potential flow has no drag.

### XFOIL solver

```bash
cd xfoil_viscous_optimizer
python run.py
```

Same prompts and same table as the in-house solver, plus the **objective**: `1` = minimum Cd at a target Cl (asks the target Cl), `2` = maximum Cl with a Cd limit (asks the maximum Cd, default 0.02). In the Cl column, `Failed` means XFOIL did not converge at the target angle, and `Timeout` means it took more than 30 s. Both are treated as failed evaluations.

### Validation

```bash
cd validation_inhouse_vs_xfoil
python run_validation.py
```

It runs NACA 0012, 2412 and 4412 over an angle-of-attack sweep (default −4° to 10°, step 2°) and explains **why** the in-house Cl differs from XFOIL.

The reference is XFOIL viscous at the real Mach number. The in-house solver differs from it for three reasons, so for every angle XFOIL is run four times (inviscid and viscous, each at Mach 0 and at the real Mach) and the difference is split into three parts:

```
dCl total = in-house Cl − XFOIL viscous Cl (real Mach)
          = num + Mach + visc
```

| Part | What it measures |
|------|------------------|
| `num` | Numerical error of the panel method: in-house vs XFOIL inviscid at Mach 0 (same physics) |
| `Mach` | Compressibility, which the in-house solver ignores |
| `visc` | Boundary layer, which the in-house solver ignores |

A positive value means the in-house solver overestimates Cl. Mach and viscosity affect each other, so their parts are the average of the two possible orders (Mach first, viscosity first); the three parts always add up exactly to the total.

If one of the intermediate XFOIL runs fails, the script uses the runs that worked instead of dropping the point, and marks the line:

| Failed run | Breakdown | Mark |
|------------|-----------|------|
| viscous at Mach 0 | Mach first | `(Mach first)` |
| inviscid at the real Mach | viscosity first | `(visc first)` |
| both of the above | Mach and viscosity together | `(Mach+visc combined)` |
| inviscid at Mach 0 | only the total | `(breakdown N/A)` |

With a single order, the split between Mach and viscosity can differ from the average by a few thousandths of Cl. In the error plot these points have hollow markers.

Example output:

```
  > Alpha =  4.0 deg ... [OK] In-House Cl: 0.7359 | XFOIL Cl: 0.6892 | dCl total +0.0467 = num -0.0078 | Mach -0.0098 | visc +0.0643
```

Typically the boundary layer explains most of the difference: it makes the real airfoil behave as if it had less camber.

---

## Output

Each optimization run creates its own folder. The `Results/` folders are created by the scripts and ignored by git, so your runs never end up in a commit.

```
Results/Results_Re<Reynolds>_Alpha<angle>_Cl<target>/       (target-Cl runs)
Results/Results_Re<Reynolds>_Alpha<angle>_CdMax<limit>/    (XFOIL, maximum-Cl runs)
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
| `export/airfoil_NACA_xxxx.dxf` | CAD: closed 2D profile in mm, ready to extrude |
| `export/airfoil_NACA_xxxx_mm.csv` | CAD: list of points `x,y,z` in mm (z = 0), trailing edge left open |
| `export/airfoil_NACA_xxxx_surface.vtk` | ParaView: airfoil contour in m with Cp and V/V∞ (in-house panel method) |
| `export/airfoil_NACA_xxxx_wing.stl` | OpenFOAM, CAD, 3D printing: the airfoil extruded along z (closed surface, in m) |

The validation script saves in `Results/Validation_Re<Reynolds>_Mach<Mach>/`:

| File | Content |
|------|---------|
| `validation_results.csv` | All Cl values (in-house and the four XFOIL runs), viscous Cd, the error parts and the breakdown method used |
| `validation_plot_cl.svg` | Cl vs alpha for each airfoil: in-house, XFOIL inviscid (Mach 0), XFOIL viscous |
| `validation_plot_error.svg` | Error parts vs alpha for each airfoil, with the total |

---

## Using the exported files

The files in `export/` use the chord you entered, so they already have the real size.

**CAD**
- The easiest route is the **DXF**: import or insert it into a sketch on the plane you want. It is in millimetres and already closed (a spline through the points plus a short line at the trailing edge), so you can extrude it straight away. For a tapered or twisted wing, place scaled or rotated copies on parallel planes and use a loft.
- The **CSV** is for CAD tools that build a curve from a list of points. It is comma-separated, in millimetres, without a header. The trailing edge is left open: close it with a line in the sketch. Check which separator and units your CAD expects before importing.
- The **STL** opens as a mesh body: fine for 3D printing or as a reference, but not editable like a solid.

**ParaView**
- Open the `.vtk` file and press Apply. Colour by `Cp` or `V_over_Vinf`.
- For a Cp–x/c plot use the filter *Plot Data* with `Points_X` on the x axis.
- The values come from the in-house panel method (potential flow), also in XFOIL runs.

**OpenFOAM**
- Copy the `.stl` into `constant/triSurface/` of your case and use it as geometry in `snappyHexMeshDict`. Coordinates are in metres, the wing goes from z = −span/2 to +span/2.
- For a 2D case, make the background mesh thinner in z than the span, so that the wing crosses the whole domain.
- Tested with OpenFOAM v1912: `surfaceCheck` reports a closed surface, one zone, consistent normals and no illegal triangles; `snappyHexMesh` removes the cells inside the airfoil and ends with "Finished meshing without any errors". In our quick test the case was not set up as a proper 2D case, so `checkMesh` then complains about the empty patches: that is the test case, not the geometry.

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
- **The in-house solver cannot minimise drag yet.** Only the XFOIL solver has the minimum-Cd objective; the in-house one needs a boundary-layer model first.
- **Many shapes give the same Cl.** The in-house search only matches the target Cl (plus the bounding box), so the airfoil it returns is one valid answer among many. The global search is random, so two runs can return different airfoils unless you enter the same seed.
- **The NACA name is rounded.** The optimizer works with continuous values (for example `m = 0.0190`, `p = 0.529`, `t = 0.1975` is saved as "NACA 2520"). The saved geometry uses the exact values, which are printed at the end of the run and written in the `.dat` header and in the CSV.
- **XFOIL does not always converge.** Airfoils where it fails get a large penalty and the search moves on.
- **Only NACA 4-digit airfoils.**
- **The export does not mesh or run a CFD case.** It gives geometry and surface data; mesh, boundary conditions and solver setup are up to you.

## Roadmap

- Boundary-layer model for the in-house solver, to estimate drag and transition without XFOIL

---

## Reference

This section explains everything you see on screen and every number the code uses. When a value comes from the original code and was never tuned, we say so.

### Symbols

| Symbol | Meaning |
|--------|---------|
| `m` | Maximum camber, as a fraction of the chord |
| `p` | Position of the maximum camber, as a fraction of the chord |
| `t` | Maximum thickness, as a fraction of the chord |
| NACA `mptt` | Name of the airfoil: first digit = m·100, second = p·10, last two = t·100, each rounded. The name is only a label: the geometry uses the exact values |
| Alpha (α) | Angle of attack, in degrees |
| Cl | Lift coefficient |
| Cd | Drag coefficient. 1 drag count = 0.0001 of Cd |
| Cp | Pressure coefficient. In the panel method Cp = 1 − (V/V∞)² |
| ΔCp | Cp,lower − Cp,upper, both taken at the same x/c. It shows where the lift is produced along the chord |
| x/c | Position along the chord: 0 = leading edge, 1 = trailing edge |
| Re | Reynolds number = speed · chord / kinematic viscosity |
| Mach | Speed / speed of sound |
| Ncrit | XFOIL transition parameter (eⁿ method): the higher it is, the later the boundary layer becomes turbulent |
| BB | Bounding box: the maximum total height the airfoil may have, in metres. It is the distance between the highest and the lowest point of the airfoil (thickness plus camber), measured at 0° angle of attack. Typical use: the space available inside a wing. `OUT` in the BB column means the airfoil is taller than this limit |
| Seed | Number that fixes the random part of the search. Same seed and same inputs give the same run |

### Terminal messages

| Prefix | Meaning |
|--------|---------|
| `[+]` | A step starts or ends |
| `[i]` | Information (computed conditions, exact parameters, seed, where files are saved) |
| `[!]` | Warning or error: read it |
| `PHASE 1 / 2 / 3` | Main steps: 1 = optimization, 2 = geometry and plots, 3 = pressure analysis of the final airfoil |
| `--> Phase 1 / Phase 2` | The two parts of the search inside PHASE 1: Genetic Algorithm, then SLSQP |

### Optimization table

```
| Eval |   m    |   p    |   t    |   Cl    |   Cd    |  BB  |    Score    |
```

| Column | Meaning |
|--------|---------|
| Eval | Number of the evaluation. The final message calls them "iterations": it is the total number of airfoils tried, including the ones outside the box or failed |
| m, p, t | Parameters of the airfoil being tried |
| Cl, Cd | Result of the analysis. Cd is `-` in the in-house solver (potential flow has no drag) |
| BB | `OUT` = the airfoil is taller than the bounding box. It is not analysed and gets a penalty |
| Score | The number the optimizer tries to make as small as possible (see below). Lower is better |

Special values in the Cl column:

| Value | Meaning |
|-------|---------|
| `Failed` | The analysis did not give a valid result. For XFOIL this includes "did not converge at the target angle" |
| `Timeout` | XFOIL took more than 30 s and was stopped |
| `-` | Not analysed (airfoil outside the bounding box) |

### Score

| Solver / objective | Score | Unit |
|--------------------|-------|------|
| In-house | (10 · (Cl − Cl_target))² | none |
| XFOIL, minimum Cd | Cd · 10⁴ + (10⁴ · Cl excess)², where Cl excess = how far Cl is outside `target ± 0.005` | drag counts |
| XFOIL, maximum Cl | −Cl + (10⁴ · Cd excess)², where Cd excess = how far Cd is above the limit | none (negative is normal) |
| Failed or timeout | 10⁹ + a small term that is lower for airfoils closer to t = 0.12 and m = 0.05 | none |
| Outside the box | 10⁹ + ((height − box) · 1000)² | none |

### Validation line

```
  > Alpha =  4.0 deg ... [OK] In-House Cl: 0.7359 | XFOIL Cl: 0.6892 | dCl total +0.0467 = num -0.0078 | Mach -0.0098 | visc +0.0643
```

| Part | Meaning |
|------|---------|
| `[OK]` | All five analyses worked |
| `[FAILED <run>]`, `[TIMEOUT <run>]` | That XFOIL run failed. `<run>` is `inviscid M0`, `inviscid M`, `viscous M0` or `viscous M` (M0 = Mach 0, M = real Mach) |
| XFOIL Cl | The reference: XFOIL viscous at the real Mach |
| dCl total | In-house Cl − reference. Positive = the in-house solver overestimates Cl |
| num | Numerical error of the panel method (in-house vs XFOIL inviscid, both at Mach 0) |
| Mach | Part due to compressibility, which the in-house solver ignores |
| visc | Part due to the boundary layer, which the in-house solver ignores |
| `(Mach first)`, `(visc first)` | A single order was used because an intermediate run failed |
| `(Mach+visc combined)` | Mach and viscosity could only be computed together |
| `(breakdown N/A)` | Only the total is available |

### Values and why

| Value | Where | Why |
|-------|-------|-----|
| Air: ν = 1.46·10⁻⁵ m²/s, a = 340.3 m/s | Fluid choice 1 | Standard atmosphere at sea level (15 °C) |
| Water: ν = 1.00·10⁻⁶ m²/s, a = 1482 m/s | Fluid choice 2 | Water at about 20 °C |
| 50 m/s, 1 m, 4°, Cl 0.8 | Default inputs | Example values from the original code. In air they give Re ≈ 3.4·10⁶ and Mach ≈ 0.15 |
| Bounding box 0.3 m | Default input | From the original code. With a 1 m chord the tallest airfoil in the search range is 0.27 m high, so by default the box never cuts anything |
| 160 panels | Default input | Accurate enough (panel-method error about 1% of Cl against XFOIL inviscid) and fast |
| 0 ≤ m ≤ 0.09 | Search range | Keeps the first digit of the NACA name a single digit |
| 0.1 ≤ p ≤ 0.7, 0.05 ≤ t ≤ 0.25 | Search range | From the original code. They cover the usual NACA 4-digit airfoils; not tuned |
| NACA 2412 | Starting airfoil, always in the initial population | Common reference airfoil, from the original code |
| Population 5, 5 generations | Genetic Algorithm | From the original code: 15 airfoils × 6 generations = 90 evaluations. A short global search is enough because SLSQP refines the result |
| SLSQP: max 50 iterations, ftol 10⁻⁴ | Local search | From the original code; not tuned |
| Step 10⁻⁴ (in-house), 2·10⁻³ (XFOIL) | SLSQP finite differences | XFOIL writes Cl with 4 decimals. With a step of 10⁻⁴ the change of Cl due to p and t is below 0.0001, so the gradient came out as zero and the search stopped early. The panel method has no such limit |
| Weight 10 | In-house score | From the original code. It only scales the score: the best airfoil does not change |
| 10⁹ | Score of failed and out-of-box airfoils | Must be higher than any valid airfoil. The worst valid score is about 4·10⁸ (Cl 2.0 away from the target), so 10⁹ always ranks a valid airfoil first |
| Small term towards t = 0.12, m = 0.05 | Failed airfoils | From the original code: among failed airfoils it prefers usual shapes, so the search moves back towards shapes that converge |
| (excess · 1000)² | Out-of-box penalty | Grows with the excess height, so the search knows which way to go back |
| ± 0.005 | Cl tolerance, minimum-Cd objective | 0.6% of a Cl of 0.8: well above XFOIL's resolution (0.0001) and small enough to keep the target meaningful. An exact Cl would be impossible to hit numerically |
| 10⁴ | Cl penalty weight, minimum-Cd objective | 0.001 of Cl outside the tolerance costs 100 drag counts. Near the target, candidate airfoils differ by about 20 counts in our runs, so going outside the tolerance never pays off |
| 10⁴ | Cd penalty weight, maximum-Cl objective | 1 drag count above the limit costs as much as 1.0 of Cl |
| Cd 0.02 | Default limit, maximum-Cl objective | From the old "maximum tolerable Cd" prompt. At 4° it is too loose to matter (the best airfoil had Cd ≈ 0.007): still to be decided |
| Ncrit = 9 | XFOIL | XFOIL's default, the usual value for an average wind tunnel |
| ITER 500 | XFOIL | Maximum viscous iterations per angle, from the original code |
| 1° steps up to the target angle | XFOIL | XFOIL converges more easily when each angle starts from the previous solution. From the original code |
| 0.1° | Angle tolerance | If XFOIL's closest converged angle is further than this from the target, its Cl belongs to another angle and the evaluation counts as failed |
| 30 s | XFOIL timeout | From the original code. A normal XFOIL run takes about a second |
| 1 – 999 999 | Random seed | Any integer works; the range just keeps it short to type |
| 10⁻⁹ | Baseline tolerance (`tests/`) | Far above round-off (we see about 10⁻¹⁵) and far below a real change (a deliberate small error in the code changed Cp by 2·10⁻⁵) |
| NACA 0012, 2412, 4412 | Validation airfoils | Symmetric, mild camber, higher camber |
| mm (DXF, CSV), m (VTK, STL) | Export units | CAD programs work in millimetres; OpenFOAM and ParaView in SI units |
| Span = chord | Default for the STL | A square planform is a neutral starting point; it only sets the length of the extrusion |
| 399 points | DXF, CSV, STL | Same points as the `.dat` file. The DXF spline stays within 0.01 mm of the exact NACA shape for a 250 mm chord |
| Straight line / flat faces at the trailing edge | DXF, STL | Closes the small gap of the NACA formula so the profile can be extruded and the STL is watertight. The geometry itself is not changed |
| −4° to 10°, step 2° | Validation sweep | From the original code |

### Conditions

- **No valid result:** if no airfoil was analysed successfully, the run stops with `[!] No airfoil could be analysed successfully` instead of reporting a meaningless airfoil.
- **Final checks (XFOIL):** after the search, XFOIL runs once more on the best airfoil. A warning appears if it does not converge, if it converges at another angle, if Cl is outside the tolerance (minimum-Cd objective) or if Cd is above the limit (maximum-Cl objective).
- **Missing libraries:** the script asks before installing them with pip. Answer no and it prints the command to run yourself.
- **Invalid input:** the optimizers stop. The validation script instead falls back to fixed values (Re 10⁶, Mach 0, alpha −4° to 10°), which are not the same as its prompt defaults.
- **Trailing edge:** the standard NACA formula leaves a small gap at the trailing edge (0.0025 of the chord for t = 0.12). The code keeps it as it is; whether to close it is still to be decided.

## Credits

- XFOIL is developed by Mark Drela and Harold Youngren (MIT) and released under the GNU General Public License.
