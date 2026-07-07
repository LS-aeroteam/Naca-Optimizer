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

def save_airfoil_coordinates(X, Y, filename):
    """
    Saves airfoil coordinates to a file in the format expected by XFOIL.

    Args:
        X (np.ndarray): X coordinates.
        Y (np.ndarray): Y coordinates.
        filename (str): The path to the output file.
    """
    with open(filename, "w") as f:
        for i in range(len(X)):
            f.write(f"{X[i]:.6f} {Y[i]:.6f}
")

