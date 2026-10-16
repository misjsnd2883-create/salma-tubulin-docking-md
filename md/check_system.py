#!/usr/bin/env python3
"""Sanity checks on a CHARMM-GUI/OpenMM system before it is simulated.

Checks, from the files in <root> (the CHARMM-GUI folder that holds `gromacs/` and `openmm/`):
  * atom count agrees between step3_input.gro and step3_input.psf (optionally against an expected value)
  * box size (optional expected values)
  * net charge of the system is zero
  * hydrogen-mass repartitioning: at least 99 % of the protein hydrogens and of the ligand hydrogens
    are heavier than 2 Da
  * time steps: 1 fs in step4_equilibration.inp and 4 fs in step5_production.inp
  * step5_production.inp contains the expected settings (310.15 K, Force-switch, PME)

Usage:
    python check_system.py --root <.../charmm-gui-XXXX> --ligand-resname LIG --expect-atoms 77701 \\
        --expect-box 124 82 82
Exit status 0 if every check passes, 1 otherwise.
"""
import argparse
import os
import sys

WATER_AND_IONS = ("TIP3", "SOD", "CLA", "POT", "MG", "CAL")


def read_inp(path):
    kv = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if "=" in line:
                key, val = line.split("=", 1)
                kv[key.strip()] = val.strip()
    return kv


def read_psf(path):
    """Atoms of a CHARMM PSF as tuples (segid, resname, atom name, charge, mass) and the NATOM header."""
    atoms, natom, inside = [], None, False
    with open(path) as fh:
        for line in fh:
            if "!NATOM" in line:
                natom, inside = int(line.split()[0]), True
                continue
            if inside:
                if not line.strip() or "!" in line:
                    break
                f = line.split()
                atoms.append((f[1], f[3], f[4], float(f[6]), float(f[7])))
    return atoms, natom


def check_system(root, expect_atoms=None, expect_box=None, ligand="LIG", log=print):
    omm = os.path.join(root, "openmm")
    ok = True

    def verdict(good, message):
        nonlocal ok
        ok &= bool(good)
        log(("PASS  " if good else "FAIL  ") + message)

    with open(os.path.join(root, "gromacs", "step3_input.gro")) as fh:
        gro = fh.read().splitlines()
    n_gro = int(gro[1])
    box = [float(x) * 10 for x in gro[-1].split()[:3]]
    atoms, n_psf = read_psf(os.path.join(omm, "step3_input.psf"))
    verdict(n_gro == n_psf == len(atoms), f"atom count: gro {n_gro:,} / psf header {n_psf:,} / psf lines {len(atoms):,}")
    if expect_atoms is not None:
        verdict(n_gro == expect_atoms, f"atom count equals the expected {expect_atoms:,}")
    log(f"      box (A): {box[0]:.1f} x {box[1]:.1f} x {box[2]:.1f}")
    if expect_box is not None:
        verdict(all(abs(b - e) < 1.5 for b, e in zip(box, expect_box)), f"box within 1.5 A of {tuple(expect_box)}")

    charge = sum(a[3] for a in atoms)
    verdict(abs(charge) < 1e-3, f"net charge {charge:+.4f}")

    solute_segs = [a for a in atoms if a[0].startswith("PRO")]
    groups = {
        "protein": [a for a in solute_segs if a[2].startswith("H")] if solute_segs else
                   [a for a in atoms if a[2].startswith("H") and a[1] not in WATER_AND_IONS + (ligand,)],
        ligand: [a for a in atoms if a[1] == ligand and a[2].startswith("H")],
    }
    for label, hs in groups.items():
        if not hs:
            verdict(False, f"{label}: no hydrogens found")
            continue
        masses = [a[4] for a in hs]
        heavy = sum(m > 2.0 for m in masses)
        verdict(heavy / len(masses) >= 0.99,
                f"{label} hydrogens: n={len(hs)}, {heavy} heavier than 2 Da, mass range {min(masses):.3f}-{max(masses):.3f}")

    for fname, want in (("step4_equilibration.inp", 0.001), ("step5_production.inp", 0.004)):
        dt = float(read_inp(os.path.join(omm, fname))["dt"])
        verdict(abs(dt - want) < 1e-12, f"{fname}: dt = {dt} ps (expected {want})")
    with open(os.path.join(omm, "step5_production.inp")) as fh:
        text = fh.read()
    for key in ("310.15", "Force-switch", "PME"):
        verdict(key in text, f"step5_production.inp contains '{key}'")
    kv = read_inp(os.path.join(omm, "step5_production.inp"))
    log("      production settings: " + ", ".join(f"{k}={kv[k]}" for k in ("nstep", "dt", "temp") if k in kv))
    return ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", required=True)
    p.add_argument("--ligand-resname", default="LIG")
    p.add_argument("--expect-atoms", type=int)
    p.add_argument("--expect-box", type=float, nargs=3)
    a = p.parse_args(argv)
    good = check_system(a.root, a.expect_atoms, a.expect_box, a.ligand_resname)
    print("\nsystem OK" if good else "\nSTOP - at least one check failed")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
