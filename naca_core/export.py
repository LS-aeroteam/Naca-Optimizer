"""
Export of the final airfoil for CAD, ParaView and OpenFOAM.

Files written in <results>/export/:
    <name>_mm.csv        x, y, z points in mm (z = 0), no header. Points only: the trailing edge is left open
    <name>.dxf           2D profile in mm, closed: spline through the points + trailing-edge line
    <name>_surface.vtk   ParaView: airfoil contour in m with Cp and V/Vinf (in-house panel method)
    <name>_wing.stl      airfoil extruded along z (m), closed binary STL meshed for surface remeshing:
                         Fluent Meshing, snappyHexMesh, CAD, ParaView

CSV, DXF and VTK use the same points as the .dat file, scaled by the chord. The STL removes section
points that are too close (cosine spacing clusters them at the trailing edge) within a set deviation.
The trailing-edge gap of the NACA formula is closed with a straight segment (DXF) or flat faces (STL).
"""
import os

import numpy as np


# ----------------------------------------------------------------------------- helpers
def _polygon_area(x, y):
    """Signed area of a closed polygon (positive = counter-clockwise)."""
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


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
# STL mesh controls, as a fraction of the chord (with chord = 1 m: 10 mm, 10 mm, 6 mm, 1.5 mm, 0.1 mm).
STL_DEFAULTS = {
    "span_element": 0.010,   # target element size along the span
    "cap_grid": 0.010,       # spacing of the interior points of the end caps (away from the outline)
    "cap_margin": 0.006,     # distance of those points from the outline (scaled down where the grading
                             # gives a finer spacing, near the leading and trailing edge)
    "min_edge": 0.0015,      # section edges shorter than this are removed (if max_deviation allows it)
    "max_deviation": 0.0001, # maximum distance of any original section point from the simplified outline
}
STL_CORNER_ANGLE = 30.0      # [deg] outline vertices that turn more than this (trailing-edge corners) are kept
STL_CAP_GRADING = 0.3        # growth of the cap point spacing with the distance from the outline


def _seg_dist(Q, A, B):
    """Distance of points Q (m, 2) from the segment A-B."""
    AB = B - A
    L2 = float(AB @ AB)
    t = np.zeros(len(Q)) if L2 == 0.0 else np.clip((Q - A) @ AB / L2, 0.0, 1.0)
    return np.linalg.norm(Q - (A + t[:, None] * AB), axis=1)


def _dist_to_outline(Q, P, chunk=2048):
    """Distance of points Q (m, 2) from the closed polyline P (n, 2)."""
    A, AB = P, np.roll(P, -1, axis=0) - P
    L2 = np.einsum("ij,ij->i", AB, AB)
    out = np.empty(len(Q))
    for s in range(0, len(Q), chunk):
        q = Q[s:s + chunk, None, :] - A[None]
        t = np.clip(np.einsum("mnj,nj->mn", q, AB) / L2, 0.0, 1.0)
        out[s:s + chunk] = np.linalg.norm(q - t[..., None] * AB, axis=2).min(axis=1)
    return out


def _span_deviation(P, a, b):
    """Max distance of the original points strictly between a and b (cyclic) from the segment P[a]-P[b]."""
    n = len(P)
    inner = [(a + k) % n for k in range(1, (b - a) % n)]
    return float(_seg_dist(P[inner], P[a], P[b]).max()) if inner else 0.0


def _simplify_section(P, min_edge, max_deviation):
    """
    Removes outline points that make edges shorter than min_edge, one at a time (shortest edge first).
    A point is removed only if every original point it represented stays within max_deviation of
    the new outline. Sharp corners (trailing edge) are never removed.
    Returns the kept indices of P (in order) and the maximum deviation.
    """
    d_in = P - np.roll(P, 1, axis=0)
    d_out = np.roll(P, -1, axis=0) - P
    turn = np.degrees(np.abs(np.arctan2(d_in[:, 0] * d_out[:, 1] - d_in[:, 1] * d_out[:, 0],
                                        np.einsum("ij,ij->i", d_in, d_out))))
    protected = turn > STL_CORNER_ANGLE

    keep = list(range(len(P)))
    while len(keep) > 4:
        Q = P[keep]
        m = len(keep)
        edges = np.linalg.norm(np.roll(Q, -1, axis=0) - Q, axis=1)
        choice = None
        for e in np.argsort(edges):
            if edges[e] >= min_edge:
                break
            candidates = []
            for k in (e, (e + 1) % m):
                if not protected[keep[k]]:
                    dev = _span_deviation(P, keep[k - 1], keep[(k + 1) % m])
                    if dev <= max_deviation:
                        candidates.append((dev, k))
            if candidates:
                choice = min(candidates)[1]
                break
        if choice is None:
            break
        del keep[choice]

    max_dev = max(_span_deviation(P, keep[i], keep[(i + 1) % len(keep)]) for i in range(len(keep)))
    return keep, max_dev


def _lattice(lo, hi, s):
    """Points of an equilateral triangular lattice with spacing s covering the box lo-hi."""
    rows = []
    for r, y in enumerate(np.arange(lo[1], hi[1] + s, s * np.sqrt(3.0) / 2.0)):
        x = np.arange(lo[0] + (r % 2) * s / 2.0, hi[0] + s, s)
        rows.append(np.column_stack([x, np.full(len(x), y)]))
    return np.vstack(rows)


def _cap_points(B, path, grid, margin):
    """
    Interior points of the end cap, graded from the outline spacing to `grid`.
    Target size at p: f(p) = min(grid, min_i(h_i + STL_CAP_GRADING * |p - B_i|)), h_i = outline spacing at B_i.
    Levels s = grid, grid/2, grid/4, ... each give a triangular lattice; level s fills where s <= f < 2s
    (the coarsest level wherever f = grid), at least margin/grid * s from the outline and at least
    0.75 s from the points already accepted (finer levels and outline). With a uniform outline spacing
    equal to grid this is a plain lattice at `grid` spacing, `margin` from the outline.
    """
    from scipy.spatial import cKDTree

    e = np.linalg.norm(np.roll(B, -1, axis=0) - B, axis=1)
    h = 0.5 * (e + np.roll(e, 1))
    tree_b = cKDTree(B)
    k = min(16, len(B))
    levels = [grid]
    while levels[-1] / 2.0 >= 0.75 * h.min():
        levels.append(levels[-1] / 2.0)
    lo, hi = B.min(axis=0), B.max(axis=0)
    accepted = [B]
    for s in reversed(levels):                       # fine to coarse
        C = _lattice(lo, hi, s)
        C = C[path.contains_points(C)]
        if len(C) == 0:
            continue
        d, idx = tree_b.query(C, k=k)
        d, idx = d.reshape(len(C), -1), idx.reshape(len(C), -1)
        f = np.minimum(grid, (h[idx] + STL_CAP_GRADING * d).min(axis=1))
        C = C[(f >= s) & ((f < 2.0 * s) | (s == grid))]
        if len(C) == 0:
            continue
        C = C[_dist_to_outline(C, B) >= margin / grid * s]
        if len(C) == 0:
            continue
        near = cKDTree(np.vstack(accepted)).query(C, k=1)[0]
        C = C[near >= 0.75 * s]
        if len(C):
            accepted.append(C)
    return np.vstack(accepted[1:]) if len(accepted) > 1 else np.empty((0, 2))


def _cap_triangulation(B, grid, margin, max_passes=30):
    """
    End-cap triangulation: outline points B (counter-clockwise) + graded interior points (_cap_points),
    Delaunay, triangles whose centroid is inside the section.
    Delaunay does not always keep the outline edges (e.g. where the airfoil is thinner than an edge):
    every missing outline edge is split at its midpoint (no geometric change) and the step is repeated.
    Returns (B with any inserted midpoints, interior points G, counter-clockwise triangles on [B; G]).
    """
    from matplotlib.path import Path
    from scipy.spatial import Delaunay

    path = Path(np.vstack([B, B[:1]]))
    G = _cap_points(B, path, grid, margin)

    for _ in range(max_passes):
        nb = len(B)
        pts = np.vstack([B, G])
        tri = Delaunay(pts).simplices
        tri = tri[path.contains_points(pts[tri].mean(axis=1))]
        # counter-clockwise orientation of every triangle
        a, b, c = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
        cross = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
        tri[cross < 0] = tri[cross < 0][:, [0, 2, 1]]
        # each outline edge must belong to exactly one cap triangle
        e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
        keys, counts = np.unique(e[:, 0] * len(pts) + e[:, 1], return_counts=True)
        count = dict(zip(keys.tolist(), counts.tolist()))
        i = np.arange(nb)
        j = (i + 1) % nb
        ob = np.minimum(i, j) * len(pts) + np.maximum(i, j)
        missing = [k for k in range(nb) if count.get(int(ob[k]), 0) != 1]
        if not missing:
            cap_area = 0.5 * np.abs(cross).sum()
            outline_area = _polygon_area(B[:, 0], B[:, 1])
            if abs(cap_area - outline_area) > 1e-9 * outline_area:
                raise RuntimeError("STL end cap: triangles do not cover the section exactly")
            return B, G, tri
        mids = 0.5 * (B[missing] + B[(np.array(missing) + 1) % nb])
        B = np.insert(B, np.array(missing) + 1, mids, axis=0)
    raise RuntimeError("STL end cap: outline edges not recovered by the triangulation")


def _mesh_quality(V, F):
    """Checks a closed indexed triangle mesh and returns its metrics. Raises RuntimeError on a defect."""
    from scipy.spatial import cKDTree

    T = V[F]
    n = len(V)
    # each edge in exactly 2 faces, with opposite directions (closed and consistently oriented)
    d = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    _, counts = np.unique(np.sort(d, axis=1) @ np.array([n, 1]), return_counts=True)
    directed_unique = len(np.unique(d @ np.array([n, 1]))) == len(d)
    if not (np.all(counts == 2) and directed_unique):
        raise RuntimeError("STL: the surface is not closed or not consistently oriented")
    if len(cKDTree(V).query_pairs(1e-9)) or len(np.unique(F.ravel())) != n:
        raise RuntimeError("STL: duplicate or unused vertices")
    cr = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    area = 0.5 * np.linalg.norm(cr, axis=1)
    if area.min() <= 0.0:
        raise RuntimeError("STL: zero-area triangle")
    volume = np.einsum("ij,ij->i", T[:, 0], np.cross(T[:, 1], T[:, 2])).sum() / 6.0
    if volume <= 0.0:
        raise RuntimeError("STL: normals point inwards (negative volume)")
    L = np.linalg.norm(T[:, [1, 2, 0]] - T[:, [2, 0, 1]], axis=2)       # L[:, k] opposite vertex k
    cosang = (L[:, [1, 2, 0]] ** 2 + L[:, [2, 0, 1]] ** 2 - L ** 2) / (2 * L[:, [1, 2, 0]] * L[:, [2, 0, 1]])
    min_angle = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))).min(axis=1)
    return {
        "vertices": n, "triangles": len(F), "volume": volume, "area": area.sum(),
        "edge_face_counts": dict(zip(*np.unique(counts, return_counts=True))),
        "min_angle": min_angle.min(), "median_angle": float(np.median(min_angle)),
        "below_10deg": int((min_angle < 10.0).sum()),
        "bbox": (V.min(axis=0), V.max(axis=0)),
    }


def write_stl_wing(X, Y, chord, span, path, name="airfoil", verbose=True, **mesh):
    """
    Airfoil extruded along z from -span/2 to +span/2, in metres, as a closed binary STL.

    The mesh is built for surface remeshing (Fluent Meshing, snappyHexMesh):
      - section points closer than min_edge are removed while the outline stays within max_deviation;
      - the side skin is split along the span into stations spaced about span_element;
      - the end caps are triangulated with interior points (Delaunay), not with outline points only;
      - vertices are shared, every edge is in 2 faces, normals point outwards (checked before writing).
    Mesh controls (keyword arguments, fraction of the chord): see STL_DEFAULTS.
    Returns a dictionary with the section changes and the mesh metrics.
    """
    unknown = set(mesh) - set(STL_DEFAULTS)
    if unknown:
        raise TypeError(f"Unknown STL mesh option(s): {sorted(unknown)}")
    opt = {k: v * chord for k, v in {**STL_DEFAULTS, **mesh}.items()}

    # Section in metres, counter-clockwise, without repeated points
    P = np.column_stack([np.asarray(X, float), np.asarray(Y, float)]) * chord
    P = P[np.linalg.norm(P - np.roll(P, 1, axis=0), axis=1) > 1e-12 * chord]
    if _polygon_area(P[:, 0], P[:, 1]) < 0:
        P = P[::-1]

    keep, max_dev = _simplify_section(P, opt["min_edge"], opt["max_deviation"])
    B = P[keep]
    n_short = int((np.linalg.norm(np.roll(B, -1, axis=0) - B, axis=1) < opt["min_edge"]).sum())
    B, G, cap = _cap_triangulation(B, opt["cap_grid"], opt["cap_margin"])
    nb, ng = len(B), len(G)

    # Stations along the span. Straight extrusion: every station is the same section at height z
    # (taper, twist or sweep would be a 2D transform of the points in `station`).
    n_span = max(1, int(np.ceil(span / opt["span_element"] - 1e-9)))
    zs = np.linspace(-span / 2.0, span / 2.0, n_span + 1)

    def station(pts2d, z):
        return np.column_stack([pts2d, np.full(len(pts2d), z)])

    V = np.vstack([station(B, z) for z in zs] + [station(G, zs[0]), station(G, zs[-1])])

    # Side skin: two triangles per outline segment and span interval (outwards for a CCW outline)
    k, i = np.meshgrid(np.arange(n_span), np.arange(nb), indexing="ij")
    a, b = k * nb + i, k * nb + (i + 1) % nb
    c, d = b + nb, a + nb
    side = np.vstack([np.stack([a, b, c], -1).reshape(-1, 3), np.stack([a, c, d], -1).reshape(-1, 3)])

    # End caps: outline indices -> first/last ring, interior indices -> interior points of that cap
    first_int = (n_span + 1) * nb

    def cap_index(ring, interior_offset):
        return np.where(cap < nb, ring * nb + cap, first_int + interior_offset + cap - nb)

    top = cap_index(n_span, ng)                 # z = +span/2, faces +z
    bottom = cap_index(0, 0)[:, [0, 2, 1]]      # z = -span/2, faces -z
    F = np.vstack([side, top, bottom])

    V32 = V.astype(np.float32)                  # values written to the file
    stats = _mesh_quality(V32.astype(np.float64), F)
    stats.update(section_points=len(P), removed_points=len(P) - len(keep), max_deviation=max_dev,
                 inserted_points=nb - len(keep), short_edges_left=n_short, interior_points=ng,
                 span_stations=n_span + 1, outline=B)

    T = V32[F]
    nrm = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
    rec = np.zeros(len(F), dtype=[("n", "<f4", (3,)), ("v", "<f4", (3, 3)), ("attr", "<u2")])
    rec["n"], rec["v"] = nrm, T
    with open(path, "wb") as f:
        f.write(f"{name} (m)".encode("ascii", "replace")[:80].ljust(80, b" "))
        f.write(np.uint32(len(F)).tobytes())
        rec.tofile(f)

    if verbose:
        print(f"      STL section: {stats['removed_points']} of {len(P)} points removed "
              f"(edges < {opt['min_edge'] * 1e3:.3g} mm), max deviation {max_dev * 1e3:.4f} mm")
        if n_short:
            print(f"      [!] {n_short} section edge(s) still shorter than {opt['min_edge'] * 1e3:.3g} mm "
                  f"(sharp corners or deviation limit)")
        print(f"      STL mesh: {stats['vertices']} vertices, {stats['triangles']} triangles, "
              f"{n_span} span elements, min angle {stats['min_angle']:.1f} deg")
    return stats

# ----------------------------------------------------------------------------- all
def export_all(X, Y, panel_results, chord, span, results_dir, name, stl_options=None):
    """
    Writes the four export files in <results_dir>/export/ and returns their paths.
    stl_options: optional STL mesh controls (fraction of the chord), see STL_DEFAULTS.
    """
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
    write_stl_wing(X, Y, chord, span, paths["stl"], name=name, **(stl_options or {}))
    return paths
