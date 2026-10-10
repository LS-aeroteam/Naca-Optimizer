"""
Check of the STL export (naca_core/export.py: write_stl_wing).

Run from anywhere:
    python tests/stl_export_check.py                 # STL files in Results/stl_export_check/
    python tests/stl_export_check.py --out <folder>

The reference case is the NACA 2421 of the run Re3424658_Alpha4.0_Cl0.8_inhouse_seed657086
(exact m, p, t of that run, 200 points per side as in run.py, chord 1 m, span 1 m).
Its metrics must be inside the thresholds below. The other cases are only reported.
Each STL is read back from the file, so the metrics describe what Fluent receives.

Exit code: 0 = reference case inside every threshold, 1 = at least one threshold missed.
"""
import argparse
import os
import sys

import numpy as np

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)
sys.path.insert(0, REPO_ROOT)

from naca_core.airfoil import naca4_airfoil                    # noqa: E402
from naca_core.export import STL_DEFAULTS, write_stl_wing      # noqa: E402

NUM_POINTS = 200                       # same as run.py (export of the optimized airfoil)
REFERENCE = ("NACA_2421", 0.02122180868323805, 0.3889605016966338, 0.20613912262357414, 1.0, 1.0)
OTHER_CASES = [                        # name, m, p, t, chord [m], span [m]
    ("NACA_0012", 0.00, 0.0, 0.12, 1.0, 1.0),
    ("NACA_9510", 0.09, 0.5, 0.10, 1.0, 1.0),
    ("NACA_9125_m09_p01_t25", 0.09, 0.1, 0.25, 1.0, 1.0),
    ("NACA_9705_m09_p07_t05", 0.09, 0.7, 0.05, 1.0, 1.0),
    ("NACA_0005", 0.00, 0.0, 0.05, 1.0, 1.0),
    ("NACA_2421_c0.25_s0.6", *REFERENCE[1:4], 0.25, 0.6),
]
# Thresholds (lengths as a fraction of the chord)
MIN_ANGLE = 8.0                        # [deg]
MEDIAN_ANGLE = (25.0, 35.0)            # [deg]
MAX_DEVIATION = STL_DEFAULTS["max_deviation"]
REL_TOL = 1e-4                         # volume and area: 0.01 %


def read_binary_stl(path):
    """Returns merged vertices (float64) and faces of a binary STL (vertices merged by exact value)."""
    with open(path, "rb") as f:
        f.read(80)
        n = int(np.frombuffer(f.read(4), "<u4")[0])
        rec = np.frombuffer(f.read(), dtype=[("n", "<f4", (3,)), ("v", "<f4", (3, 3)), ("a", "<u2")], count=n)
    pts = rec["v"].reshape(-1, 3)
    V, inv = np.unique(pts, axis=0, return_inverse=True)
    return V.astype(np.float64), inv.reshape(-1, 3), rec["n"].astype(np.float64)


def outline_at_top(V, F):
    """Ordered outline of the top cap: edges at z = max shared with a side face."""
    ztop = V[:, 2].max()
    at_top = np.isclose(V[:, 2], ztop)
    side = F[~np.all(at_top[F], axis=1)]
    nxt = {}
    for a, b, c in side:
        for p, q in ((a, b), (b, c), (c, a)):
            if at_top[p] and at_top[q]:
                nxt[p] = q
    start = next(iter(nxt))
    loop = [start]
    while nxt[loop[-1]] != start:
        loop.append(nxt[loop[-1]])
        if len(loop) > len(nxt):
            raise RuntimeError("top outline is not a single loop")
    return V[loop, :2]


def dist_to_outline(Q, P):
    A, AB = P, np.roll(P, -1, axis=0) - P
    L2 = np.einsum("ij,ij->i", AB, AB)
    q = Q[:, None, :] - A[None]
    t = np.clip(np.einsum("mnj,nj->mn", q, AB) / L2, 0.0, 1.0)
    return np.linalg.norm(q - t[..., None] * AB, axis=2).min(axis=1)


def metrics(path, X, Y, chord, span):
    from scipy.spatial import cKDTree

    V, F, N = read_binary_stl(path)
    T = V[F]
    n = len(V)
    d = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    _, counts = np.unique(np.sort(d, axis=1) @ np.array([n, 1]), return_counts=True)
    cr = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    area_t = 0.5 * np.linalg.norm(cr, axis=1)
    L = np.linalg.norm(T[:, [1, 2, 0]] - T[:, [2, 0, 1]], axis=2)
    cosang = (L[:, [1, 2, 0]] ** 2 + L[:, [2, 0, 1]] ** 2 - L ** 2) / (2 * L[:, [1, 2, 0]] * L[:, [2, 0, 1]])
    ang = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))).min(axis=1)
    is_cap = np.ptp(T[:, :, 2], axis=1) == 0.0

    # unmodified section (the points given to the export), scaled by the chord
    P0 = np.column_stack([X, Y]) * chord
    A0 = abs(0.5 * np.sum(P0[:, 0] * np.roll(P0[:, 1], -1) - np.roll(P0[:, 0], -1) * P0[:, 1]))
    L0 = np.linalg.norm(np.roll(P0, -1, axis=0) - P0, axis=1).sum()
    bbox0 = np.array([[P0[:, 0].min(), P0[:, 1].min(), -span / 2], [P0[:, 0].max(), P0[:, 1].max(), span / 2]])
    outline = outline_at_top(V, F)
    e_out = np.linalg.norm(np.roll(outline, -1, axis=0) - outline, axis=1)

    return {
        "vertices": n,
        "triangles": len(F),
        "edges_per_face_count": {int(k): int(v) for k, v in zip(*np.unique(counts, return_counts=True))},
        "oriented": len(np.unique(d @ np.array([n, 1]))) == len(d),
        "duplicate_vertices": len(cKDTree(V).query_pairs(1e-9 * chord)),
        "zero_area_triangles": int((area_t <= 0.0).sum()),
        "normals_match_winding": bool(np.all(np.einsum("ij,ij->i", N, cr) > 0)),
        "volume": np.einsum("ij,ij->i", T[:, 0], np.cross(T[:, 1], T[:, 2])).sum() / 6.0,
        "volume_ref": A0 * span,
        "area": area_t.sum(),
        "area_ref": L0 * span + 2.0 * A0,
        "bbox": np.array([V.min(axis=0), V.max(axis=0)]),
        "bbox_ref": bbox0,
        "deviation": dist_to_outline(P0, outline).max(),
        "outline_points": len(outline),
        "outline_edge_min": e_out.min(),
        "min_angle": ang.min(),
        "min_angle_side": ang[~is_cap].min(),
        "min_angle_caps": ang[is_cap].min(),
        "median_angle": float(np.median(ang)),
        "below_10deg": int((ang < 10.0).sum()),
    }


def check(m, chord):
    """List of (label, value text, passed)."""
    rel_v = m["volume"] / m["volume_ref"] - 1.0
    rel_a = m["area"] / m["area_ref"] - 1.0
    bbox_err = np.abs(m["bbox"] - m["bbox_ref"]).max()
    closed = set(m["edges_per_face_count"]) == {2}
    return [
        ("edges: faces per edge (count)", str(m["edges_per_face_count"]), closed),
        ("edges: consistent orientation", str(m["oriented"]), m["oriented"]),
        ("duplicate vertices", str(m["duplicate_vertices"]), m["duplicate_vertices"] == 0),
        ("zero-area triangles", str(m["zero_area_triangles"]), m["zero_area_triangles"] == 0),
        ("normals consistent with winding", str(m["normals_match_winding"]), m["normals_match_winding"]),
        ("volume [m^3] (ref, diff)", f"{m['volume']:.7f} ({m['volume_ref']:.7f}, {rel_v * 100:+.5f} %)",
         m["volume"] > 0 and abs(rel_v) <= REL_TOL),
        ("area [m^2] (ref, diff)", f"{m['area']:.6f} ({m['area_ref']:.6f}, {rel_a * 100:+.5f} %)",
         abs(rel_a) <= REL_TOL),
        ("bbox min [m]", np.array2string(m["bbox"][0], precision=6), True),
        ("bbox max [m]", np.array2string(m["bbox"][1], precision=6), True),
        ("bbox max diff from section [mm]", f"{bbox_err * 1e3:.4f}", bbox_err <= MAX_DEVIATION * chord),
        ("max deviation [mm]", f"{m['deviation'] * 1e3:.4f}", m["deviation"] <= MAX_DEVIATION * chord),
        ("min angle [deg] (side, caps)",
         f"{m['min_angle']:.2f} ({m['min_angle_side']:.2f}, {m['min_angle_caps']:.2f})",
         m["min_angle"] >= MIN_ANGLE),
        ("median of min angle [deg]", f"{m['median_angle']:.2f}",
         MEDIAN_ANGLE[0] <= m["median_angle"] <= MEDIAN_ANGLE[1]),
        ("triangles < 10 deg", f"{m['below_10deg']} of {m['triangles']}", True),
    ]


def run_case(case, out_dir):
    name, mp, pp, tp, chord, span = case
    X, Y, _ = naca4_airfoil(mp, pp, tp, num_points=NUM_POINTS)
    path = os.path.join(out_dir, f"airfoil_{name}_wing.stl")
    print(f"\n=== {name}  (m={mp:.4f} p={pp:.4f} t={tp:.4f}, chord {chord} m, span {span} m)")
    write_stl_wing(X, Y, chord, span, path, name=f"airfoil_{name}")
    m = metrics(path, X, Y, chord, span)
    print(f"  file: {path}")
    print(f"  {m['vertices']} vertices, {m['triangles']} triangles, outline {m['outline_points']} points "
          f"(shortest edge {m['outline_edge_min'] * 1e3:.3f} mm)")
    rows = check(m, chord)
    for label, text, ok in rows:
        print(f"  {'ok  ' if ok else 'FAIL'} {label:34s} {text}")
    return all(ok for _, _, ok in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default=os.path.join(REPO_ROOT, "Results", "stl_export_check"))
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("STL mesh controls (fraction of chord):", STL_DEFAULTS)
    print(f"Thresholds: min angle >= {MIN_ANGLE} deg, median {MEDIAN_ANGLE[0]}-{MEDIAN_ANGLE[1]} deg, "
          f"deviation <= {MAX_DEVIATION} c, volume/area within {REL_TOL * 100} %")
    print("\n##### Reference case (thresholds enforced)")
    passed = run_case(REFERENCE, args.out)
    print("\n##### Other cases (report only)")
    others = {case[0]: run_case(case, args.out) for case in OTHER_CASES}

    print("\n##### Summary")
    print(f"  {REFERENCE[0]:28s} {'PASS' if passed else 'FAIL'} (enforced)")
    for name, ok in others.items():
        print(f"  {name:28s} {'inside thresholds' if ok else 'outside thresholds'} (report only)")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
