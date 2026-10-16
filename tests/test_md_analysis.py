"""Synthetic-data tests for analysis/md_analysis.py.

Each test builds a small protein-ligand system with known geometry, so the expected
numbers are known analytically. Run:  python -m pytest tests/   (or python tests/test_md_analysis.py)
"""
import os
import sys
import tempfile

import mdtraj as md
import numpy as np
from mdtraj.core import element as E
from scipy.spatial.transform import Rotation

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import md_analysis as A  # noqa: E402

LIG_ATOMS = [  # name, element, position (nm) relative to the ligand origin
    ("C1", E.carbon, (0.00, 0.00, 0.00)), ("C2", E.carbon, (0.15, 0.00, 0.00)),
    ("C3", E.carbon, (0.30, 0.00, 0.00)), ("N1", E.nitrogen, (0.45, 0.00, 0.00)),
    ("O1", E.oxygen, (0.60, 0.00, 0.00)), ("O2", E.oxygen, (0.75, 0.00, 0.00)),
    ("F1", E.fluorine, (0.90, 0.00, 0.00)), ("C4", E.carbon, (0.30, 0.15, 0.00)),
    ("H1", E.hydrogen, (0.00, 0.00, 0.25)), ("H2", E.hydrogen, (0.15, 0.00, 0.25)),
    ("H3", E.hydrogen, (0.30, 0.00, 0.25)), ("H4", E.hydrogen, (0.45, 0.00, 0.25)),
]
PROT_ATOMS = [("N", E.nitrogen), ("CA", E.carbon), ("C", E.carbon), ("O", E.oxygen), ("CB", E.carbon)]
BOX_NM = 60.0


class System:
    """Two protein chains of ALA-like residues parked on a 3 nm grid, plus one ligand."""

    def __init__(self, chain_sizes=(32, 32), lig="LIG", names=None):
        names = names or {}
        self.top = md.Topology()
        self.res_atoms = []
        n = 0
        for ci, size in enumerate(chain_sizes):
            chain = self.top.add_chain()
            for k in range(size):
                nm, seq = names.get(n, ("ALA", k + 2))
                r = self.top.add_residue(nm, chain, resSeq=seq)
                self.res_atoms.append({a: self.top.add_atom(a, el, r).index for a, el in PROT_ATOMS})
                n += 1
        lchain = self.top.add_chain()
        lres = self.top.add_residue(lig, lchain, resSeq=1)
        self.lig = {a: self.top.add_atom(a, el, lres).index for a, el, _ in LIG_ATOMS}
        self.n_atoms = self.top.n_atoms
        self.n_res = n
        side = int(np.ceil(n ** (1 / 3)))
        grid = [(i, j, k) for i in range(side) for j in range(side) for k in range(side)]
        self.park = [np.array(grid[i], float) * 3.0 + 10.0 for i in range(n)]  # nm, far from the ligand

    def base_xyz(self):
        xyz = np.zeros((self.n_atoms, 3))
        for r, atoms in enumerate(self.res_atoms):
            for a, idx in atoms.items():
                xyz[idx] = self.park[r]
        for a, _, pos in LIG_ATOMS:
            xyz[self.lig[a]] = np.array(pos) + 1.0  # ligand origin at (1,1,1) nm
        return xyz

    def traj(self, xyz_frames):
        t = md.Trajectory(np.asarray(xyz_frames, dtype=np.float32), self.top)
        t.unitcell_lengths = np.full((t.n_frames, 3), BOX_NM, dtype=np.float32)
        t.unitcell_angles = np.full((t.n_frames, 3), 90.0, dtype=np.float32)
        return t


def lig_origin():
    return np.array([1.0, 1.0, 1.0])


def lig_pos(name):
    return lig_origin() + np.array([p for n, _, p in LIG_ATOMS if n == name][0])


# ------------------------------------------------------------------ tests
def test_rmsd_alignment_and_heavy_atoms_only():
    s = System()
    n = 100
    frames = np.repeat(s.base_xyz()[None], n, axis=0)
    # ligand heavy atoms shift 0.1 nm (1 A) along x from frame 50; hydrogens shift an extra 0.5 nm
    heavy = [i for a, i in s.lig.items() if not a.startswith("H")]
    hyd = [i for a, i in s.lig.items() if a.startswith("H")]
    frames[50:, heavy, 0] += 0.1
    frames[50:, hyd, 0] += 0.6
    # random rigid-body motion of the whole system in every frame (rotation + translation)
    rng = np.random.default_rng(1)
    for f in range(n):
        R = Rotation.random(random_state=int(rng.integers(1 << 30))).as_matrix()
        frames[f] = frames[f] @ R.T + rng.normal(0, 2.0, 3)
    lig_rmsd, prot_rmsd, n_heavy = A.rmsd_series(s.traj(frames), "LIG")
    assert n_heavy == 8, n_heavy
    assert np.allclose(lig_rmsd[:50], 0.0, atol=0.01), lig_rmsd[:50].max()
    assert np.allclose(lig_rmsd[50:], 1.0, atol=0.01), (lig_rmsd[50:].min(), lig_rmsd[50:].max())
    assert prot_rmsd.max() < 0.01, prot_rmsd.max()  # protein does not move internally


def test_rmsf_matches_sigma_sqrt3():
    s = System()
    n = 2000
    rng = np.random.default_rng(2)
    sigma, noisy_res, sigma_big = 0.05, 7, 0.15  # nm per coordinate
    frames = np.repeat(s.base_xyz()[None], n, axis=0)
    for r, atoms in enumerate(s.res_atoms):
        frames[:, atoms["CA"], :] += rng.normal(0, sigma_big if r == noisy_res else sigma, (n, 3))
    for f in range(n):  # plus rigid-body motion, which the alignment has to remove
        R = Rotation.random(random_state=int(rng.integers(1 << 30))).as_matrix()
        frames[f] = frames[f] @ R.T + rng.normal(0, 2.0, 3)
    labels, rmsf = A.ca_rmsf(s.traj(frames))
    assert len(rmsf) == s.n_res and len(labels) == s.n_res
    expected = np.full(s.n_res, sigma * np.sqrt(3) * 10)  # Angstrom
    expected[noisy_res] = sigma_big * np.sqrt(3) * 10
    rel = np.abs(rmsf - expected) / expected
    assert rel.max() < 0.06, rel.max()
    assert np.argmax(rmsf) == noisy_res


def test_contact_occupancy_known_geometry():
    names = {5: ("ALA", 6), 20: ("SER", 20), 30: ("GLY", 30), 33: ("LYS", 352), 40: ("ASP", 360), 50: ("VAL", 99)}
    s = System(names=names)
    n = 100
    xyz = s.base_xyz()
    frames = np.repeat(xyz[None], n, axis=0)
    O = lambda r: s.res_atoms[r]["O"]
    frames[:, O(5)] = lig_pos("C1") + (-0.30, 0, 0)                      # always in contact           -> 1.0
    frames[:50, O(20)] = lig_pos("C1") + (0, 0.36, 0)                    # contact in first half       -> 0.5
    frames[50:, O(20)] = lig_pos("C1") + (0, 0.46, 0)
    frames[:, O(30)] = lig_pos("C1") + (0, 0, 0.55)                      # only the H atom is near     -> 0.0
    frames[:, O(33)] = lig_pos("N1") + (0, -0.39, 0)                     # 3.9 A, inside the cutoff    -> 1.0
    frames[:, O(40)] = lig_pos("F1") + (0, 0, -0.41)                     # 4.1 A, outside the cutoff   -> 0.0
    occ = A.contact_occupancy(s.traj(frames), "LIG")
    assert abs(occ["αAla6"] - 1.0) < 1e-9, occ
    assert abs(occ["αSer20"] - 0.5) < 1e-9, occ
    assert abs(occ["βLys352"] - 1.0) < 1e-9, occ
    for absent in ("αGly30", "βAsp360", "βVal99"):
        assert occ.get(absent, 0.0) == 0.0, (absent, occ.get(absent))


def test_polar_contacts_known_geometry():
    names = {33: ("LYS", 352), 34: ("ASP", 353), 35: ("ALA", 354), 36: ("GLY", 355), 37: ("SER", 356)}
    s = System(names=names)
    n = 100
    frames = np.repeat(s.base_xyz()[None], n, axis=0)
    at = lambda r, a: s.res_atoms[r][a]
    frames[:40, at(33, "N")] = lig_pos("O1") + (0, 0.34, 0)              # 3.4 A for 40 % of frames -> 0.4
    frames[40:, at(33, "N")] = lig_pos("O1") + (0, 0.60, 0)
    frames[:, at(34, "O")] = lig_pos("O2") + (0, 0, 0.36)                # 3.6 A, outside cutoff     -> absent
    frames[:, at(35, "CB")] = lig_pos("O1") + (0, -0.30, 0)              # carbon, not a polar atom  -> absent
    frames[:, at(36, "O")] = lig_pos("C1") + (0, -0.30, 0)               # near a ligand carbon only -> absent
    frames[:, at(37, "N")] = lig_pos("N1") + (0, 0.30, 0)                # 3.0 A to ligand N         -> 1.0
    occ = A.polar_contact_occupancy(s.traj(frames), "LIG")
    assert set(occ) == {"βLys352:N", "βSer356:N"}, occ
    assert abs(occ["βLys352:N"] - 0.4) < 1e-9
    assert abs(occ["βSer356:N"] - 1.0) < 1e-9


def test_labels_and_chains():
    s = System(names={33: ("LYS", 352)})
    t = s.traj(np.repeat(s.base_xyz()[None], 2, axis=0))
    labels = A.chain_labels(t)
    assert list(labels.values()) == ["α", "β"]
    res = [r for r in t.topology.residues if r.name == "LYS"][0]
    assert A.residue_label(res, labels) == "βLys352"


def _write_rep(rep, s, n_seg, frames_per_seg, encode=True):
    os.makedirs(rep, exist_ok=True)
    xyz = s.base_xyz()
    md.Trajectory(xyz[None].astype(np.float32), s.top).save_pdb(os.path.join(rep, "step3_input.pdb"))
    open(os.path.join(rep, "step5_production.inp"), "w").write("nstep = 250000  # steps\ndt = 0.004  # ps\n")
    for i in range(1, n_seg + 2):  # one extra, unfinished segment without a .done marker
        fr = np.repeat(xyz[None], frames_per_seg, axis=0)
        if encode:
            fr[:, s.res_atoms[0]["CA"], 0] += i * 0.001  # segment number encoded in the coordinates
        s.traj(fr).save_dcd(os.path.join(rep, f"step5_{i}.dcd"))
        if i <= n_seg:
            open(os.path.join(rep, f"step5_{i}.done"), "w").write("ok")


def test_segment_order_time_axis_and_unfinished_segment_ignored():
    s = System()
    with tempfile.TemporaryDirectory() as d:
        _write_rep(d, s, n_seg=12, frames_per_seg=2)  # 12 segments: step5_10 must come after step5_9
        assert A.completed_segments(d) == list(range(1, 13))
        traj, t, seg_ns, n_seg = A.load_production(d, "LIG")
        assert seg_ns == 1.0 and n_seg == 12 and traj.n_frames == 24
        assert np.isclose(t[0], 0.5) and np.isclose(t[-1], 12.0)
        ca = s.res_atoms[0]["CA"]
        keep = traj.topology.select("protein or resname LIG")
        col = list(keep).index(ca)
        code = np.round((traj.xyz[:, col, 0] - s.base_xyz()[ca, 0]) / 0.001).astype(int)
        assert list(code[::2]) == list(range(1, 13)), code  # frames are in segment order


def test_cli_end_to_end():
    s = System()
    with tempfile.TemporaryDirectory() as d:
        _write_rep(d, s, n_seg=3, frames_per_seg=10, encode=False)
        out = os.path.join(d, "out.npz")
        A.main(["--rep-dir", d, "--ligand-resname", "LIG", "--out", out])
        z = np.load(out)
        for key in ("time_ns", "ligand_rmsd", "protein_rmsd", "rmsf", "rmsf_labels",
                    "contact_labels", "contact_occ", "polar_labels", "polar_occ"):
            assert key in z.files, key
        assert len(z["time_ns"]) == 30 and np.isclose(z["time_ns"][-1], 3.0)
        assert z["ligand_rmsd"].max() < 0.01


def test_plot_smoke():
    import plot_md as P
    s = System()
    with tempfile.TemporaryDirectory() as d:
        rep = os.path.join(d, "rep")
        _write_rep(rep, s, n_seg=3, frames_per_seg=10, encode=False)
        out = os.path.join(d, "x.npz")
        A.main(["--rep-dir", rep, "--ligand-resname", "LIG", "--out", out])
        P.main(["--s12", out, "--cn2", out, "--outdir", os.path.join(d, "fig")])
        for f in ("rmsd_compare.png", "rmsf_compare.png"):
            assert os.path.getsize(os.path.join(d, "fig", f)) > 10_000, f


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print("PASS", t.__name__)
    print(f"{len(tests)} tests passed")
