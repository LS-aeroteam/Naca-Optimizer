import numpy as np

def naca4_airfoil(m_param, p_param, t_param, chord=1.0, num_points=100):
    """
    Generates the coordinates for a NACA 4-digit airfoil.
    Uses cosine spacing for XFOIL compatibility.

    Args:
        m_param (float): Maximum camber as a fraction of the chord (e.g., 0.02 for NACA 2412).
        p_param (float): Position of maximum camber as a fraction of the chord (e.g., 0.4 for NACA 2412).
        t_param (float): Maximum thickness as a fraction of the chord (e.g., 0.12 for NACA 2412).
        chord (float, optional): The airfoil chord length. Defaults to 1.0.
        num_points (int, optional): The number of points for the upper or lower surface. Defaults to 100.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: X coordinates of the airfoil surface, starting from the trailing edge,
                          going along the upper surface to the leading edge, and back along the lower surface.
            - np.ndarray: Y coordinates of the airfoil surface.
            - tuple: A detailed breakdown of (xu, yu, xl, yl) coordinates.
    """
    beta = np.linspace(0, np.pi, num_points)
    x = chord * (0.5 * (1 - np.cos(beta)))
    
    # Thickness distribution
    yt = 5 * t_param * chord * (
        0.2969 * np.sqrt(x / chord)
        - 0.1260 * (x / chord)
        - 0.3516 * (x / chord) ** 2
        + 0.2843 * (x / chord) ** 3
        - 0.1015 * (x / chord) ** 4
    )

    if p_param == 0 or m_param == 0:  # Symmetric airfoil
        xu, yu = x, yt
        xl, yl = x, -yt
    else:  # Cambered airfoil
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
            
        theta = np.arctan(dyc_dx)
        xu = x - yt * np.sin(theta)
        yu = yc + yt * np.cos(theta)
        xl = x + yt * np.sin(theta)
        yl = yc - yt * np.cos(theta)

    # Combine to a single closed contour
    X = np.concatenate((np.flip(xu), xl[1:]))
    Y = np.concatenate((np.flip(yu), yl[1:]))
    
    return X, Y, (xu, yu, xl, yl)

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

    def __init__(self, key, params, descriptions, bounds, initial_guess, coordinates, name, failure_score):
        """
        Args:
            coordinates (callable): (params, num_points) -> (X, Y), closed contour TE -> upper -> LE -> lower -> TE.
            name (callable): params -> airfoil name (e.g. "NACA 2412").
            failure_score (callable): (params, penalty) -> score of a failed analysis: the penalty plus a small
                term, lower for usual shapes, so that the search moves back towards shapes that work.
        """
        self.key = key
        self.params = tuple(params)
        self.descriptions = tuple(descriptions)
        self.bounds = list(bounds)
        self.initial_guess = list(initial_guess)
        self._coordinates = coordinates
        self._name = name
        self._failure_score = failure_score

    def coordinates(self, params, num_points):
        """Airfoil contour (X, Y) with num_points points per surface."""
        return self._coordinates(params, num_points)

    def name(self, params):
        return self._name(params)

    def failure_score(self, params, penalty):
        return self._failure_score(params, penalty)

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
    coordinates=lambda params, num_points: naca4_airfoil(*params, num_points=num_points)[:2],
    name=_naca4_name,
    # From the original code: among failed airfoils it prefers t = 0.12 and m = 0.05
    failure_score=lambda params, penalty: penalty + (0.12 - params[2]) ** 2 * 1e5 + (0.05 - params[0]) ** 2 * 1e5,
)


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