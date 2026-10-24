"""
Run inputs: defaults, fluids and the Python input file (python run.py --input inputs.py).

Only the standard library is used here: run.py imports this module before
perform_all_checks() has verified numpy/scipy/matplotlib.

The input file is a plain Python file with one UPPERCASE variable per input (SPEED = 50.0, ...).
load_input_file() returns the same dictionary as the interactive questions in run.py,
plus the solver, so the rest of run.py does not depend on where the inputs come from.
"""
import runpy

INHOUSE = "inhouse"
XFOIL = "xfoil"

# Objectives, as written in the input file. "1"/"2" are the answers of the interactive question.
OBJECTIVE_CHOICES = {"min_cd": "1", "max_cl": "2"}

# Airfoil families (naca_core.airfoil.get_family). Chosen only in the input file:
# the interactive questions always use "naca4".
AIRFOIL_FAMILIES = ("naca4", "naca5")

FLUIDS = {
    "air": {"name": "Air (Standard SL)", "viscosity": 1.46e-5, "speed_of_sound": 340.3},
    "water": {"name": "Water (20 degrees)", "viscosity": 1.00e-6, "speed_of_sound": 1482.0},
}

# Default of every input: used by the interactive questions (Enter) and for the variables
# missing from the input file (written there in UPPERCASE). None = random seed / span = chord.
DEFAULTS = {
    "solver": INHOUSE,
    "airfoil_family": "naca4",
    "fluid": "air",
    "speed": 50,
    "chord": 1.0,
    "alpha": 4.0,
    "objective": "min_cd",
    "target_cl": 0.8,
    "max_cd": 0.006,
    "max_height_box": 0.3,
    "num_panels": 160,
    "seed": None,
    "span": None,
}

# Expected type of every input. seed and span can also be None.
_NUMBER = "a number"
_INTEGER = "a whole number"
_TEXT = "text in quotes"
_KEY_TYPES = {
    "solver": _TEXT, "airfoil_family": _TEXT, "fluid": _TEXT, "speed": _NUMBER, "chord": _NUMBER,
    "alpha": _NUMBER, "objective": _TEXT, "target_cl": _NUMBER, "max_cd": _NUMBER, "max_height_box": _NUMBER,
    "num_panels": _INTEGER, "seed": _INTEGER, "span": _NUMBER,
}
_CAN_BE_NONE = ("seed", "span")
_TEXT_CHOICES = {"solver": (INHOUSE, XFOIL), "airfoil_family": AIRFOIL_FAMILIES, "fluid": tuple(FLUIDS),
                 "objective": tuple(OBJECTIVE_CHOICES)}


class InputFileError(Exception):
    """Invalid input file: the message says which variable and why."""


def flow_conditions(fluid, speed, chord):
    """Reynolds and Mach number from the fluid properties, the speed and the chord."""
    reynolds = (speed * chord) / fluid["viscosity"]
    mach = speed / fluid["speed_of_sound"]
    return reynolds, mach


def _check_type(key, value):
    name = key.upper()
    if value is None and key in _CAN_BE_NONE:
        return
    expected = _KEY_TYPES[key]
    # bool is a subclass of int in Python: "CHORD = True" must not be read as 1
    if isinstance(value, bool):
        ok = False
    elif expected == _NUMBER:
        ok = isinstance(value, (int, float))
    elif expected == _INTEGER:
        ok = isinstance(value, int)
    else:
        ok = isinstance(value, str)
    if not ok:
        allowed = f"{expected} or None" if key in _CAN_BE_NONE else expected
        raise InputFileError(f"{name} must be {allowed}, got {value!r}")
    if key in _TEXT_CHOICES and value not in _TEXT_CHOICES[key]:
        choices = " | ".join(f'"{c}"' for c in _TEXT_CHOICES[key])
        raise InputFileError(f"{name} must be {choices}, got {value!r}")


def load_input_file(path, solver=None):
    """
    Runs a Python input file and reads its UPPERCASE variables. Every variable is optional:
    a missing one takes its value from DEFAULTS.

    Args:
        path (str): Path of the input file.
        solver (str, optional): Solver chosen on the command line; it replaces the one in the file.

    Returns:
        dict: solver and airfoil_family, plus the same keys as run.get_user_input (reynolds, mach,
              alpha, objective_choice, target_cl, max_cd, max_height_box, chord, num_panels, seed, span).
              seed is None when the file does not set it (run.py then draws a random one).

    The values are not checked beyond what the interactive questions check: numbers must be
    numbers, and solver, fluid and objective must be one of their choices.

    Raises:
        InputFileError: file missing, Python error in the file, unknown variable, wrong type or choice.
    """
    try:
        namespace = runpy.run_path(path)
    except FileNotFoundError:
        raise InputFileError("file not found") from None
    except SyntaxError as e:
        raise InputFileError(f"Python syntax error at line {e.lineno}: {e.msg}") from None
    except Exception as e:
        raise InputFileError(f"error while running the file: {type(e).__name__}: {e}") from None

    # Only UPPERCASE names are inputs: lowercase helpers and imports in the file are ignored
    data = {name.lower(): value for name, value in namespace.items()
            if name.isupper() and not name.startswith("_")}
    unknown = sorted(name.upper() for name in set(data) - set(DEFAULTS))
    if unknown:
        valid = ", ".join(key.upper() for key in DEFAULTS)
        raise InputFileError(f"unknown variable(s): {', '.join(unknown)} (valid names: {valid})")
    for key, value in data.items():
        _check_type(key, value)

    values = {**DEFAULTS, **data}
    if solver is not None:
        values["solver"] = solver
    reynolds, mach = flow_conditions(FLUIDS[values["fluid"]], values["speed"], values["chord"])
    # Same rule as the interactive questions: the objective is chosen only with XFOIL,
    # the in-house solver (no Cd) always reaches the target Cl
    objective_choice = OBJECTIVE_CHOICES[values["objective"]] if values["solver"] == XFOIL else "1"
    return {
        'solver': values["solver"], 'airfoil_family': values["airfoil_family"],
        'reynolds': reynolds, 'mach': mach, 'alpha': float(values["alpha"]),
        'objective_choice': objective_choice,
        'target_cl': float(values["target_cl"]) if objective_choice == "1" else None,
        'max_cd': float(values["max_cd"]) if objective_choice == "2" else None,
        'max_height_box': float(values["max_height_box"]), 'chord': float(values["chord"]),
        'num_panels': values["num_panels"], 'seed': values["seed"],
        'span': float(values["span"]) if values["span"] is not None else float(values["chord"]),
    }
