"""
Export of the final airfoil for CAD, ParaView and OpenFOAM.

Files written in <results>/export/:
    <name>_mm.csv        x, y, z points in mm (z = 0), no header. Points only: the trailing edge is left open
    <name>.dxf           2D profile in mm, closed: spline through the points + trailing-edge line
    <name>_surface.vtk   ParaView: airfoil contour in m with Cp and V/Vinf (in-house panel method)
    <name>_wing.stl      airfoil extruded along z (m), closed surface: CAD, ParaView, OpenFOAM snappyHexMesh

Only the file format changes: the points are the same as in the .dat file, scaled by the chord.
The trailing-edge gap of the NACA formula is closed with a straight segment (DXF) or flat faces (STL).
"""
import os

import numpy as np


# ----------------------------------------------------------------------------- helpers
def _polygon_area(x, y):
    """Signed area of a closed polygon (positive = counter-clockwise)."""
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def _triangulate_polygon(x, y):
    """
    Ear-clipping triangulation of a simple polygon (no holes).

    Returns a list of (i, j, k) index triples, all counter-clockwise.
    Needed for the STL end caps: the airfoil outline is not convex when the camber is high,
    so a simple fan from one point would give overlapping triangles.
    """
    n = len(x)
    idx = list(range(n))
    if _polygon_area(x, y) < 0:
        idx.reverse()

    def cross(a, b, c):
        return (x[b] - x[a]) * (y[c] - y[a]) - (y[b] - y[a]) * (x[c] - x[a])

    def inside(p, a, b, c):
        # point p strictly inside triangle abc (counter-clockwise)
        return cross(a, b, p) > 0 and cross(b, c, p) > 0 and cross(c, a, p) > 0

    triangles = []
    guard = 0
    while len(idx) > 3:
        ear_found = False
        m = len(idx)
        for k in range(m):
            a, b, c = idx[k - 1], idx[k], idx[(k + 1) % m]
            if cross(a, b, c) <= 0:
                continue  # reflex or collinear vertex: not an ear
            if any(inside(p, a, b, c) for p in idx if p not in (a, b, c)):
                continue
            triangles.append((a, b, c))
            del idx[k]
            ear_found = True
            break
        if not ear_found:
            # Only collinear vertices left: remove one and continue
            guard += 1
            if guard > n:
                raise RuntimeError("Polygon triangulation failed")
            del idx[0]
    triangles.append(tuple(idx))
    return triangles


# ----------------------------------------------------------------------------- CAD
def write_csv_mm(X, Y, chord, path):
    """x, y, z = 0 in mm, one point per line, no header (e.g. SolidWorks 'Curve Through XYZ Points')."""
    scale = chord * 1000.0
    with open(path, "w") as f:
        for x, y in zip(X, Y):
            f.write(f"{x * scale:.6f},{y * scale:.6f},0.000000\n")


def write_dxf_mm(X, Y, chord, path):
    """
    Closed 2D profile in mm: one spline through all points (control-vertex form, readable by
    most CAD programs) plus a straight line that closes the trailing edge.
    """
    import ezdxf

    scale = chord * 1000.0
    pts = [(float(x) * scale, float(y) * scale) for x, y in zip(X, Y)]

    doc = ezdxf.new("R2010", setup=False, units=4)  # units=4: millimetres
    msp = doc.modelspace()
    doc.layers.add("AIRFOIL")
    msp.add_cad_spline_control_frame(pts, dxfattribs={"layer": "AIRFOIL"})
    msp.add_line(pts[-1], pts[0], dxfattribs={"layer": "AIRFOIL"})  # trailing-edge closure
    doc.saveas(path)


# ----------------------------------------------------------------------------- ParaView
def write_vtk_surface(panel_results, chord, path, title="airfoil"):
    """
    Legacy VTK polydata (ASCII): the panel control points joined by lines, in metres,
    with Cp and V/Vinf as point data. Opens directly in ParaView.
    """
    XC = np.asarray(panel_results["XC"]) * chord
    YC = np.asarray(panel_results["YC"]) * chord
    Cp = np.asarray(panel_results["Cp"])
    v_ratio = np.sqrt(np.maximum(1.0 - Cp, 0.0))  # Cp = 1 - (V/Vinf)^2
    n = len(XC)
    n_half = panel_results["num_panels"] // 2
    surface = np.array([0] * n_half + [1] * (n - n_half))  # 0 = lower, 1 = upper

    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"{title} - in-house panel method (potential flow)\n")
        f.write("ASCII\nDATASET POLYDATA\n")
        f.write(f"POINTS {n} double\n")
        for x, y in zip(XC, YC):
            f.write(f"{x:.9e} {y:.9e} 0.0\n")
        f.write(f"LINES 1 {n + 1}\n")
        f.write(str(n) + " " + " ".join(str(i) for i in range(n)) + "\n")
        f.write(f"POINT_DATA {n}\n")
        f.write("SCALARS Cp double 1\nLOOKUP_TABLE default\n")
        f.write("\n".join(f"{v:.9e}" for v in Cp) + "\n")
        f.write("SCALARS V_over_Vinf double 1\nLOOKUP_TABLE default\n")
        f.write("\n".join(f"{v:.9e}" for v in v_ratio) + "\n")
        f.write("SCALARS surface_0lower_1upper int 1\nLOOKUP_TABLE default\n")
        f.write("\n".join(str(v) for v in surface) + "\n")


# ----------------------------------------------------------------------------- OpenFOAM / 3D
def write_stl_wing(X, Y, chord, span, path, name="airfoil"):
    """
    Airfoil extruded along z from -span/2 to +span/2, in metres, as a closed ASCII STL
    (side skin + trailing-edge strip + two end caps). Normals point outwards.
    """
    x = np.asarray(X, dtype=float) * chord
    y = np.asarray(Y, dtype=float) * chord
    n = len(x)
    z0, z1 = -span / 2.0, span / 2.0

    ccw = _polygon_area(x, y) > 0
    tris = []

    # Side skin: one quad (two triangles) per segment, trailing-edge closure included (i = n-1 -> 0)
    for i in range(n):
        j = (i + 1) % n
        a0, b0 = (x[i], y[i], z0), (x[j], y[j], z0)
        a1, b1 = (x[i], y[i], z1), (x[j], y[j], z1)
        if ccw:
            tris += [(a0, b0, b1), (a0, b1, a1)]
        else:
            tris += [(a0, b1, b0), (a0, a1, b1)]

    # End caps: z1 faces +z, z0 faces -z
    for (i, j, k) in _triangulate_polygon(x, y):
        tris.append(((x[i], y[i], z1), (x[j], y[j], z1), (x[k], y[k], z1)))
        tris.append(((x[i], y[i], z0), (x[k], y[k], z0), (x[j], y[j], z0)))

    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for p0, p1, p2 in tris:
            nrm = np.cross(np.subtract(p1, p0), np.subtract(p2, p0))
            length = np.linalg.norm(nrm)
            nrm = nrm / length if length > 0 else nrm
            f.write(f"  facet normal {nrm[0]:.6e} {nrm[1]:.6e} {nrm[2]:.6e}\n    outer loop\n")
            for p in (p0, p1, p2):
                f.write(f"      vertex {p[0]:.9e} {p[1]:.9e} {p[2]:.9e}\n")
            f.write("    endloop\n  endfacet\n")
        f.write(f"endsolid {name}\n")


# ----------------------------------------------------------------------------- all
def export_all(X, Y, panel_results, chord, span, results_dir, name):
    """Writes the four export files in <results_dir>/export/ and returns their paths."""
    export_dir = os.path.join(results_dir, "export")
    os.makedirs(export_dir, exist_ok=True)
    paths = {
        "csv": os.path.join(export_dir, f"{name}_mm.csv"),
        "dxf": os.path.join(export_dir, f"{name}.dxf"),
        "vtk": os.path.join(export_dir, f"{name}_surface.vtk"),
        "stl": os.path.join(export_dir, f"{name}_wing.stl"),
    }
    write_csv_mm(X, Y, chord, paths["csv"])
    write_dxf_mm(X, Y, chord, paths["dxf"])
    write_vtk_surface(panel_results, chord, paths["vtk"], title=name)
    write_stl_wing(X, Y, chord, span, paths["stl"], name=name)
    return paths
