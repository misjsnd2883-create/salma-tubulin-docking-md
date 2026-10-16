"""Trajectory analysis for the Salma-12 and CN2 (control) MD runs.

Computes, from the production DCD segments of one replicate:
  * ligand heavy-atom RMSD and protein C-alpha RMSD (vs. the first production frame,
    after superposition on all protein C-alpha atoms)
  * per-residue heavy-atom contact occupancy (any ligand heavy atom within 4.0 A)
  * polar-contact occupancy (ligand N/O within 3.5 A of protein N/O, per protein atom)
  * per-residue C-alpha RMSF (about the mean structure, after superposition)

Usage:
    python md_analysis.py --rep-dir <.../rep1> --ligand-resname LIG --out salma12.npz
"""
import argparse
import glob
import os
import re
from collections import Counter

import mdtraj as md
import numpy as np

CONTACT_CUTOFF_NM = 0.40
POLAR_CUTOFF_NM = 0.35
GREEK = ("α", "β")  # alpha, beta (chain order in the CHARMM-GUI systems)


# ---------------------------------------------------------------- input files
def read_inp(path):
    """Parse 'key = value  # comment' lines of a CHARMM-GUI/OpenMM .inp file."""
    kv = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if "=" in line:
                key, val = line.split("=", 1)
                kv[key.strip()] = val.strip()
    return kv


def segment_length_ns(inp_path):
    """Length of one production segment in ns: nstep * dt (dt in ps)."""
    kv = read_inp(inp_path)
    try:
        return int(kv["nstep"]) * float(kv["dt"]) / 1000.0
    except KeyError as err:
        raise KeyError(f"{err} not found in {inp_path}; pass --segment-ns explicitly") from None


def completed_segments(rep_dir):
    """Segment numbers with a .done marker, sorted numerically (step5_10 after step5_9)."""
    ids = []
    for f in glob.glob(os.path.join(rep_dir, "step5_*.done")):
        m = re.fullmatch(r"step5_(\d+)\.done", os.path.basename(f))
        if m:
            ids.append(int(m.group(1)))
    return sorted(ids)


def is_hydrogen(atom):
    sym = atom.element.symbol if atom.element is not None else ""
    return sym == "H" or atom.name.upper().startswith("H")


def load_production(rep_dir, ligand_resname, topology="step3_input.pdb", segment_ns=None):
    """Load protein + ligand atoms of all completed production segments.

    Returns (traj, time_ns) where time_ns[i] is the end time of frame i within the
    production run (first frame = one frame interval after the start of production).
    """
    ids = completed_segments(rep_dir)
    if not ids:
        raise FileNotFoundError(f"no completed segments (step5_*.done) in {rep_dir}")
    if segment_ns is None:
        segment_ns = segment_length_ns(os.path.join(rep_dir, "step5_production.inp"))
    top = md.load_topology(os.path.join(rep_dir, topology))
    keep = top.select(f"protein or resname {ligand_resname}")
    dcds = [os.path.join(rep_dir, f"step5_{i}.dcd") for i in ids]
    traj = md.load(dcds, top=top, atom_indices=keep)
    frames_per_segment = traj.n_frames / len(ids)
    dt_frame = segment_ns / frames_per_segment
    time_ns = dt_frame * (np.arange(traj.n_frames) + 1)
    return traj, time_ns, segment_ns, len(ids)


# --------------------------------------------------------------- atom helpers
def ligand_atoms(traj, resname, heavy_only=True):
    idx = [a.index for a in traj.topology.atoms
           if a.residue.name == resname and not (heavy_only and is_hydrogen(a))]
    if not idx:
        raise ValueError(f"no atoms with residue name {resname}")
    return np.array(idx)


def protein_heavy(traj):
    return np.array([a.index for a in traj.topology.atoms
                     if a.residue.is_protein and not is_hydrogen(a)])


def chain_labels(traj):
    """Map protein chain index -> alpha/beta symbol (first protein chain = alpha)."""
    chains = sorted({r.chain.index for r in traj.topology.residues if r.is_protein})
    return {c: (GREEK[i] if i < len(GREEK) else f"c{c}") for i, c in enumerate(chains)}


def residue_label(res, labels):
    return f"{labels[res.chain.index]}{res.name.capitalize()}{res.resSeq}"


def atom_label(atom, labels):
    return f"{residue_label(atom.residue, labels)}:{atom.name}"


# ------------------------------------------------------------------- analyses
def rmsd_series(traj, ligand_resname):
    """Ligand heavy-atom RMSD and protein C-alpha RMSD (Angstrom) vs. frame 0.

    All frames are superposed on the protein C-alpha atoms of frame 0 first; both
    RMSDs are then computed directly from the superposed coordinates (the C-alpha
    superposition is already the RMSD-optimal fit for the protein).
    """
    t = traj[:]  # copy
    ca = t.topology.select("protein and name CA")
    lig = ligand_atoms(t, ligand_resname, heavy_only=True)
    t.superpose(t, 0, atom_indices=ca)
    lig_rmsd = np.sqrt(np.mean(np.sum((t.xyz[:, lig] - t.xyz[0, lig]) ** 2, axis=2), axis=1)) * 10
    prot_rmsd = np.sqrt(np.mean(np.sum((t.xyz[:, ca] - t.xyz[0, ca]) ** 2, axis=2), axis=1)) * 10
    return lig_rmsd, prot_rmsd, len(lig)


def contact_occupancy(traj, ligand_resname, cutoff_nm=CONTACT_CUTOFF_NM):
    """Fraction of frames in which any ligand heavy atom is within the cutoff of any
    heavy atom of a protein residue. Returns {label: occupancy}."""
    lig = ligand_atoms(traj, ligand_resname, heavy_only=True)
    prot = protein_heavy(traj)
    neighbours = md.compute_neighbors(traj, cutoff_nm, lig, haystack_indices=prot)
    atoms = list(traj.topology.atoms)
    counts = Counter()
    for frame_atoms in neighbours:
        residues = {atoms[i].residue.index for i in frame_atoms}
        counts.update(residues)
    labels = chain_labels(traj)
    res = list(traj.topology.residues)
    return {residue_label(res[r], labels): n / traj.n_frames for r, n in counts.items()}


def polar_contact_occupancy(traj, ligand_resname, cutoff_nm=POLAR_CUTOFF_NM):
    """Per protein N/O atom: fraction of frames with any ligand N/O atom within the
    cutoff. Returns {label: occupancy} for atoms with occupancy > 0."""
    sym = lambda a: a.element.symbol if a.element is not None else ""
    lig = np.array([a.index for a in traj.topology.atoms
                    if a.residue.name == ligand_resname and sym(a) in ("N", "O")])
    prot = np.array([a.index for a in traj.topology.atoms
                     if a.residue.is_protein and sym(a) in ("N", "O")])
    pairs = np.array([(i, j) for i in lig for j in prot])
    d = md.compute_distances(traj, pairs).reshape(traj.n_frames, len(lig), len(prot))
    occ = (d < cutoff_nm).any(axis=1).mean(axis=0)
    labels = chain_labels(traj)
    atoms = list(traj.topology.atoms)
    return {atom_label(atoms[prot[k]], labels): float(occ[k]) for k in np.nonzero(occ)[0]}


def ca_rmsf(traj, n_iter=3):
    """Per-residue C-alpha RMSF in Angstrom, about the mean structure.

    Frames are superposed on all protein C-alpha atoms, first to frame 0, then
    iteratively to the running mean structure.
    Returns (labels, rmsf_angstrom).
    """
    ca = traj.topology.select("protein and name CA")
    t = traj.atom_slice(ca)
    t.superpose(t, 0)
    for _ in range(n_iter):
        ref = md.Trajectory(t.xyz.mean(axis=0, keepdims=True), t.topology)
        t.superpose(ref, 0)
    mean = t.xyz.mean(axis=0)
    rmsf = np.sqrt(np.mean(np.sum((t.xyz - mean) ** 2, axis=2), axis=0)) * 10
    chains = chain_labels(traj)
    labels = [f"{chains[a.residue.chain.index]}{a.residue.resSeq}" for a in t.topology.atoms]
    return labels, rmsf


# ------------------------------------------------------------------------ CLI
def top_items(d, n):
    return sorted(d.items(), key=lambda kv: -kv[1])[:n]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rep-dir", required=True)
    p.add_argument("--ligand-resname", required=True, help="LIG for Salma-12, CN2 for the control")
    p.add_argument("--out", required=True, help="output .npz")
    p.add_argument("--topology", default="step3_input.pdb")
    p.add_argument("--segment-ns", type=float, default=None,
                   help="segment length in ns (default: nstep*dt from step5_production.inp)")
    p.add_argument("--top", type=int, default=15, help="how many contacts to print")
    a = p.parse_args(argv)

    traj, time_ns, seg_ns, n_seg = load_production(a.rep_dir, a.ligand_resname, a.topology, a.segment_ns)
    chains = chain_labels(traj)
    print(f"segments: {n_seg} x {seg_ns:g} ns = {n_seg * seg_ns:g} ns | frames: {traj.n_frames} "
          f"| frame interval: {time_ns[1] - time_ns[0]:.3g} ns")
    for c, g in chains.items():
        n = sum(1 for r in traj.topology.residues if r.is_protein and r.chain.index == c)
        print(f"protein chain {c} -> {g}: {n} residues")

    lig_rmsd, prot_rmsd, n_heavy = rmsd_series(traj, a.ligand_resname)
    print(f"\nligand heavy atoms: {n_heavy}")
    print(f"ligand heavy-atom RMSD (A): mean {lig_rmsd.mean():.2f} | max {lig_rmsd.max():.2f} | last {lig_rmsd[-1]:.2f}")
    print(f"protein C-alpha RMSD  (A): mean {prot_rmsd.mean():.2f} | max {prot_rmsd.max():.2f} | last {prot_rmsd[-1]:.2f}")

    lig_idx = ligand_atoms(traj, a.ligand_resname)
    cen = traj.xyz[:, lig_idx].mean(axis=1)
    step = np.linalg.norm(np.diff(cen, axis=0), axis=1).max() * 10
    print(f"largest frame-to-frame ligand centroid displacement: {step:.2f} A"
          + ("  <-- check for periodic-boundary wrapping" if step > 10 else ""))

    contacts = contact_occupancy(traj, a.ligand_resname)
    print(f"\ncontact occupancy (heavy atoms <= {CONTACT_CUTOFF_NM * 10:.1f} A), top {a.top}:")
    for k, v in top_items(contacts, a.top):
        print(f"  {k:14s} {100 * v:5.1f}%")

    polar = polar_contact_occupancy(traj, a.ligand_resname)
    print(f"\npolar contacts (ligand N/O - protein N/O <= {POLAR_CUTOFF_NM * 10:.1f} A), top 8:")
    for k, v in top_items(polar, 8):
        print(f"  {k:20s} {100 * v:5.1f}%")

    rmsf_labels, rmsf = ca_rmsf(traj)
    print(f"\nC-alpha RMSF (A): mean {rmsf.mean():.2f} | median {np.median(rmsf):.2f} | max {rmsf.max():.2f}")

    np.savez(a.out, time_ns=time_ns, ligand_rmsd=lig_rmsd, protein_rmsd=prot_rmsd,
             rmsf_labels=np.array(rmsf_labels), rmsf=rmsf,
             contact_labels=np.array(list(contacts)), contact_occ=np.array(list(contacts.values())),
             polar_labels=np.array(list(polar)), polar_occ=np.array(list(polar.values())))
    print(f"\nsaved {a.out}")


if __name__ == "__main__":
    main()
