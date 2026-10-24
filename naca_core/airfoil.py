import numpy as np
from scipy.optimize import brentq


def _cosine_x(chord, num_points):
    """Cosine spacing: more points near the leading and trailing edge."""
    beta = np.linspace(0, np.pi, num_points)
    return chord * (0.5 * (1 - np.cos(beta)))


def _naca_thickness(x, t_param, chord, closed_te):
    """NACA 4/5-digit thickness distribution (half thickness), open or closed trailing edge."""
    # Last coefficient: -0.1015 (original, open TE) or -0.1036 (closed TE)
    a4 = 0.1036 if closed_te else 0.1015
    yt = 5 * t_param * chord * (
        0.2969 * np.sqrt(x / chord)
        - 0.1260 * (x / chord)
        - 0.3516 * (x / chord) ** 2
        + 0.2843 * (x / chord) ** 3
        - a4 * (x / chord) ** 4
    )
    if closed_te:
        # The coefficients add up to zero, so yt = 0 at x = chord: set it exactly, otherwise round-off
        # leaves about -1e-17 and the upper surface ends a hair below the lower one
        yt[-1] = 0.0
    return yt


def _contour(x, yt, yc, dyc_dx):
    """
    Thickness applied perpendicular to the mean line, then one closed contour
    TE -> upper -> LE -> lower -> TE. Returns X, Y and (xu, yu, xl, yl).
    """
    theta = np.arctan(dyc_dx)
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    X = np.concatenate((np.flip(xu), xl[1:]))
    Y = np.concatenate((np.flip(yu), yl[1:]))
    return X, Y, (xu, yu, xl, yl)


def naca4_airfoil(m_param, p_param, t_param, chord=1.0, num_points=100, closed_te=False):
    """
    Generates the coordinates for a NACA 4-digit airfoil.
    Uses cosine spacing for XFOIL compatibility.

    Trailing edge:
        open (default)  original NACA formula: a small gap at x = chord (0.0025 of the chord for
                        t = 0.12). Used for every analysis: XFOIL's viscous solution is reliable
                        on it, while on a sharp trailing edge it jumps between two solutions.
        closed          last coefficient of the thickness -0.1036 instead of -0.1015: zero thickness
                        at x = chord, upper and lower surface end at the same point. Used for the
                        exported files (CAD, mesh), which need a closed profile.

    Args:
        m_param (float): Maximum camber as a fraction of the chord (e.g., 0.02 for NACA 2412).
        p_param (float): Position of maximum camber as a fraction of the chord (e.g., 0.4 for NACA 2412).
        t_param (float): Maximum thickness as a fraction of the chord (e.g., 0.12 for NACA 2412).
        chord (float, optional): The airfoil chord length. Defaults to 1.0.
        num_points (int, optional): The number of points for the upper or lower surface. Defaults to 100.
        closed_te (bool, optional): Closed (sharp) trailing edge. Defaults to False (open).

    Returns:
        tuple: A tuple containing:
            - np.ndarray: X coordinates of the airfoil surface, starting from the trailing edge,
                          going along the upper surface to the leading edge, and back along the lower surface.
            - np.ndarray: Y coordinates of the airfoil surface.
            - tuple: A detailed breakdown of (xu, yu, xl, yl) coordinates.
    """
    x = _cosine_x(chord, num_points)
    yt = _naca_thickness(x, t_param, chord, closed_te)

    if p_param == 0 or m_param == 0:  # Symmetric airfoil
        xu, yu = x, yt
        xl, yl = x, -yt
        X = np.concatenate((np.flip(xu), xl[1:]))
        Y = np.concatenate((np.flip(yu), yl[1:]))
        return X, Y, (xu, yu, xl, yl)

    # Cambered airfoil
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)

    # Front section
    front_mask = x < p_param * chord
    if np.any(front_mask):
        x_front = x[front_mask]
        yc[front_mask] = (m_param / p_param ** 2) * (2 * p_param * (x_front / chord) - (x_front / chord) ** 2)
        dyc_dx[front_mask] = (2 * m_param / p_param ** 2) * (p_param - x_front / chord)

    # Back section
    back_mask = x >= p_param * chord
    if np.any(back_mask):
        x_back = x[back_mask]
        yc[back_mask] = (m_param / (1 - p_param) ** 2) * ((1 - 2 * p_param) + 2 * p_param * (x_back / chord) - (x_back / chord) ** 2)
        dyc_dx[back_mask] = (2 * m_param / (1 - p_param) ** 2) * (p_param - x_back / chord)

    return _contour(x, yt, yc, dyc_dx)


def naca5_mean_line_constants(cl_design, p_param):
    """
    r and k1 of the NACA 5-digit non-reflex mean line (Abbott & von Doenhoff), for any p.

    r follows from the position of maximum camber, p = r * (1 - sqrt(r / 3)); k1 from the design
    lift coefficient through the thin-airfoil relation, k1 = 6 * cl / Q(r). With the tabulated r
    they reproduce the tabulated k1 within 0.03 % for p = 0.15, 0.20, 0.25 (230, 240, 250).
    """
    r = brentq(lambda rr: rr * (1 - np.sqrt(rr / 3)) - p_param, p_param, 1.0)
    q = ((3 * r - 7 * r ** 2 + 8 * r ** 3 - 4 * r ** 4) / np.sqrt(r * (1 - r))
         - 1.5 * (1 - 2 * r) * (np.pi / 2 - np.arcsin(1 - 2 * r)))
    return r, 6 * cl_design / q


def naca5_airfoil(cl_design, p_param, t_param, chord=1.0, num_points=100, closed_te=False):
    """
    Generates the coordinates for a NACA 5-digit airfoil, non-reflex mean line (e.g. NACA 23012).

    Same thickness distribution, point spacing, trailing-edge options and output as naca4_airfoil.

    Args:
        cl_design (float): Design lift coefficient (first digit * 0.15, e.g. 0.3 for NACA 23012).
        p_param (float): Position of maximum camber as a fraction of the chord
            (second digit / 20, e.g. 0.15 for NACA 23012).
        t_param (float): Maximum thickness as a fraction of the chord (e.g. 0.12 for NACA 23012).
        chord, num_points, closed_te: as in naca4_airfoil.
    """
    x = _cosine_x(chord, num_points)
    yt = _naca_thickness(x, t_param, chord, closed_te)
    r, k1 = naca5_mean_line_constants(cl_design, p_param)
    xc = x / chord
    front = xc < r
    yc = chord * np.where(front, k1 / 6 * (xc ** 3 - 3 * r * xc ** 2 + r ** 2 * (3 - r) * xc),
                          k1 * r ** 3 / 6 * (1 - xc))
    dyc_dx = np.where(front, k1 / 6 * (3 * xc ** 2 - 6 * r * xc + r ** 2 * (3 - r)), -k1 * r ** 3 / 6)
    return _contour(x, yt, yc, dyc_dx)


CST_TE_GAP = 0.0025  # [fraction of chord] trailing-edge gap of the analysed CST airfoils, as NACA 2412


def _bernstein(x, order):
    """Bernstein polynomials of the given order at x: array (len(x), order + 1)."""
    from math import comb
    return np.column_stack([comb(order, i) * x ** i * (1 - x) ** (order - i) for i in range(order + 1)])


def _cst_class(x):
    """CST class function of a round-nose airfoil with a sharp trailing edge: sqrt(x) * (1 - x)."""
    return np.sqrt(x) * (1 - x)


def cst_airfoil(a_upper, a_lower, chord=1.0, num_points=100, te_gap=0.0):
    """
    Generates the coordinates for a CST (Kulfan) airfoil: upper and lower surface are independent,

        y(x) = sqrt(x) * (1 - x) * sum_i A_i * B_i(x)  +/-  x * te_gap / 2

    with B_i the Bernstein polynomials of order len(A) - 1. The lower-surface coefficients of a
    usual airfoil are negative. Same point spacing and output as naca4_airfoil.

    Args:
        a_upper (sequence): Upper-surface coefficients.
        a_lower (sequence): Lower-surface coefficients (same number as a_upper).
        chord (float, optional): The airfoil chord length. Defaults to 1.0.
        num_points (int, optional): The number of points for the upper or lower surface. Defaults to 100.
        te_gap (float, optional): Trailing-edge gap as a fraction of the chord. 0 = closed trailing edge.
    """
    x = _cosine_x(1.0, num_points)
    basis = _bernstein(x, len(a_upper) - 1) * _cst_class(x)[:, None]
    yu = basis @ np.asarray(a_upper, float) + x * te_gap / 2
    yl = basis @ np.asarray(a_lower, float) - x * te_gap / 2
    xu = xl = x
    X = chord * np.concatenate((np.flip(xu), xl[1:]))
    Y = chord * np.concatenate((np.flip(yu), yl[1:]))
    return X, Y, (chord * xu, chord * yu, chord * xl, chord * yl)


def cst_fit(xs, ys, order):
    """Least-squares CST coefficients of one surface given by points (xs, ys), chord 1."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    inner = (xs > 0) & (xs < 1)  # the class function is zero at both ends
    basis = _bernstein(xs[inner], order) * _cst_class(xs[inner])[:, None]
    return np.linalg.lstsq(basis, ys[inner], rcond=None)[0]


class AirfoilFamily:
    """
    A family of airfoils described by a few parameters: everything the optimizer and run.py
    need to know about the geometry, so they do not depend on one family.

    Attributes:
        key (str): Name of the family (e.g. "naca4").
        params (tuple): Parameter names, in the order of the parameter vector (e.g. ("m", "p", "t")).
        descriptions (tuple): Short description of each parameter (plot legends).
        bounds (list): Search limits (min, max) of each parameter.
        initial_guess (list): Starting airfoil, always in the initial population of the search.
    """

    def __init__(self, key, params, descriptions, bounds, initial_guess, coordinates, name, failure_score,
                 shape_error=None):
        """
        Args:
            coordinates (callable): (params, num_points, closed_te) -> (X, Y), contour TE -> upper -> LE -> lower -> TE.
            name (callable): params -> airfoil name (e.g. "NACA 2412").
            failure_score (callable): (params, penalty) -> score of a failed analysis: the penalty plus a small
                term, lower for usual shapes, so that the search moves back towards shapes that work.
            shape_error (callable, optional): params -> None for a valid shape, or a short label (shown in the
                Cl column) when the parameters give no real airfoil; the shape is then not analysed.
        """
        self.key = key
        self.params = tuple(params)
        self.descriptions = tuple(descriptions)
        self.bounds = list(bounds)
        self.initial_guess = list(initial_guess)
        self._coordinates = coordinates
        self._name = name
        self._failure_score = failure_score
        self._shape_error = shape_error

    def coordinates(self, params, num_points, closed_te=False):
        """Airfoil contour (X, Y) with num_points points per surface. closed_te: sharp trailing edge (export)."""
        return self._coordinates(params, num_points, closed_te)

    def name(self, params):
        return self._name(params)

    def failure_score(self, params, penalty):
        return self._failure_score(params, penalty)

    def shape_error(self, params):
        return self._shape_error(params) if self._shape_error else None

    def describe(self, params, digits=4):
        """Exact parameters as text, e.g. "m = 0.0200, p = 0.4000, t = 0.1200"."""
        return ", ".join(f"{n} = {v:.{digits}f}" for n, v in zip(self.params, params))


def _naca4_name(params):
    m, p, t = params
    return f"NACA {int(round(m*100))}{int(round(p*10))}{int(round(t*100)):02d}"


NACA4 = AirfoilFamily(
    key="naca4",
    params=("m", "p", "t"),
    descriptions=("Max Camber", "Camber Position", "Thickness"),
    # 0 <= m <= 0.09 keeps the first digit of the name a single digit; p and t from the original code
    bounds=[(0.0, 0.09), (0.1, 0.7), (0.05, 0.25)],
    initial_guess=[0.02, 0.4, 0.12],  # NACA 2412
    coordinates=lambda params, num_points, closed_te: naca4_airfoil(*params, num_points=num_points,
                                                                    closed_te=closed_te)[:2],
    name=_naca4_name,
    # From the original code: among failed airfoils it prefers t = 0.12 and m = 0.05
    failure_score=lambda params, penalty: penalty + (0.12 - params[2]) ** 2 * 1e5 + (0.05 - params[0]) ** 2 * 1e5,
)


def _naca5_name(params):
    cl, p, t = params
    return f"NACA {int(round(cl / 0.15))}{int(round(p * 20))}0{int(round(t * 100)):02d}"


NACA5 = AirfoilFamily(
    key="naca5",
    params=("cl_d", "p", "t"),  # cl_d: design lift coefficient (not the Cl computed by the solver)
    descriptions=("Design Cl", "Camber Position", "Thickness"),
    # p: range of the standard mean lines (210 to 250); cl up to 0.6 (first digit 4); t as NACA 4
    bounds=[(0.05, 0.6), (0.05, 0.25), (0.05, 0.25)],
    initial_guess=[0.3, 0.15, 0.12],  # NACA 23012
    coordinates=lambda params, num_points, closed_te: naca5_airfoil(*params, num_points=num_points,
                                                                    closed_te=closed_te)[:2],
    name=_naca5_name,
    # Like NACA 4: among failed airfoils it prefers usual shapes, here t = 0.12 and cl = 0.3 (NACA 23012)
    failure_score=lambda params, penalty: penalty + (0.12 - params[2]) ** 2 * 1e5 + (0.3 - params[0]) ** 2 * 1e3,
)


def _max_distance(points, polyline):
    """Largest distance of the points (n, 2) from the polyline (m, 2)."""
    A, AB = polyline[:-1], np.diff(polyline, axis=0)
    t = np.clip(np.einsum("nmj,mj->nm", points[:, None, :] - A[None], AB) / np.einsum("mj,mj->m", AB, AB), 0, 1)
    return float(np.linalg.norm(points[:, None, :] - (A[None] + t[..., None] * AB[None]), axis=2).min(axis=1).max())


def cst_family(order=3):
    """
    Free airfoil family: CST (Kulfan) with order + 1 coefficients on each surface.

    Search limits: the CST fit of NACA 4-digit airfoils spread over the whole NACA 4-digit search
    range (m 0/0.045/0.09, p 0.1/0.4/0.7, t 0.05/0.15/0.25), widened by 25 % of the spread on each
    side. Fits further than 0.005 of the chord from their NACA airfoil are left out: those NACA
    shapes (high camber close to the leading edge) reach ahead of x = 0, which a CST airfoil cannot,
    and their coefficients would only widen the limits. Starting airfoil: the CST fit of NACA 2412,
    as for NACA 4-digit. The analyses use a trailing-edge gap of CST_TE_GAP (as NACA 2412), the
    exported airfoil a closed trailing edge.
    """
    fits = []
    for m in (0.0, 0.045, 0.09):
        for p in (0.1, 0.4, 0.7):
            for t in (0.05, 0.15, 0.25):
                X, Y, (xu, yu, xl, yl) = naca4_airfoil(m, p, t, num_points=200, closed_te=True)
                coeffs = np.concatenate((cst_fit(xu, yu, order), cst_fit(xl, yl, order)))
                CX, CY, _ = cst_airfoil(coeffs[:order + 1], coeffs[order + 1:], num_points=400)
                if _max_distance(np.column_stack([X, Y])[::4], np.column_stack([CX, CY])) <= 0.005:
                    fits.append(coeffs)
    _, _, (xu, yu, xl, yl) = naca4_airfoil(0.02, 0.4, 0.12, num_points=200, closed_te=True)
    start = np.concatenate((cst_fit(xu, yu, order), cst_fit(xl, yl, order)))
    fits = np.array(fits)
    lo, hi = fits.min(axis=0), fits.max(axis=0)
    bounds = [(float(a - 0.25 * (b - a)), float(b + 0.25 * (b - a))) for a, b in zip(lo, hi)]
    n = order + 1
    spread = np.array([b - a for a, b in bounds])
    x_check = _cosine_x(1.0, 201)[1:-1]
    thickness_basis = _bernstein(x_check, order) * _cst_class(x_check)[:, None]

    def shape_error(params):
        # The upper surface must stay above the lower one (zero thickness only at the two ends)
        thickness = thickness_basis @ (np.asarray(params[:n]) - np.asarray(params[n:]))
        return "Crossed" if thickness.min() <= 0 else None

    return AirfoilFamily(
        key=f"cst{order}",
        params=tuple(f"Au{i}" for i in range(n)) + tuple(f"Al{i}" for i in range(n)),
        descriptions=tuple(f"Upper coefficient {i}" for i in range(n)) + tuple(f"Lower coefficient {i}" for i in range(n)),
        bounds=bounds,
        initial_guess=[float(v) for v in start],
        coordinates=lambda params, num_points, closed_te: cst_airfoil(
            params[:n], params[n:], num_points=num_points, te_gap=0.0 if closed_te else CST_TE_GAP)[:2],
        name=lambda params: f"CST-{order}",
        # Among failed airfoils it prefers shapes close to the starting airfoil (NACA 2412)
        failure_score=lambda params, penalty: penalty + 1e3 * float(np.sum(((np.asarray(params) - start) / spread) ** 2)),
        shape_error=shape_error,
    )


def get_family(key, cst_order=3):
    """Airfoil family from its name in the input file ("naca4", "naca5", "cst")."""
    if key == "cst":
        return cst_family(cst_order)
    return {"naca4": NACA4, "naca5": NACA5}[key]


def save_airfoil_coordinates(X, Y, filename, header=None):
    """
    Saves airfoil coordinates to a file in the format expected by XFOIL.

    Args:
        X (np.ndarray): X coordinates.
        Y (np.ndarray): Y coordinates.
        filename (str): The path to the output file.
        header (str, optional): Airfoil name written on the first line.
            XFOIL reads a non-numeric first line as the airfoil name.
    """
    with open(filename, "w") as f:
        if header:
            f.write(f"{header}\n")
        for i in range(len(X)):
            f.write(f"{X[i]:.6f} {Y[i]:.6f}\n")