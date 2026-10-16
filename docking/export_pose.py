#!/usr/bin/env python3
"""Turn a docked pose into the ligand file and the protein-ligand complex used to build an MD system.

Steps:
  1. read the pose (PDBQT written by dock.py) with Meeko; connectivity and bond orders come from
     the SMILES stored in the PDBQT file
  2. optionally deprotonate carboxylic acids (docking used the neutral acids; the simulation of
     Salma-12 used the carboxylate anion)
  3. add hydrogens with RDKit (placed from the geometry of the docked heavy atoms)
  4. write <out>.sdf (ligand with hydrogens) and <out>_complex.pdb (receptor site from
     dock.py + the ligand as residue LIG)

The CIP label of the stereocentre and the heavy-atom coordinates of the docked pose are checked
against the output, and the script stops if they differ.

Usage:
    python export_pose.py --pose work/Salma-12_s2.pdbqt --site work/1SA0_site.pdb --out work/S12 --anion
"""
import argparse
import os
import sys

from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dock import pose_to_mol  # noqa: E402

ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")


def deprotonate_acids(mol):
    """Carboxylic acid -> carboxylate (formal charge -1 on the former OH oxygen). Returns (mol, count)."""
    mol = Chem.RWMol(mol)
    hits = mol.GetSubstructMatches(ACID)
    for _c, _o, oh in hits:
        mol.GetAtomWithIdx(oh).SetFormalCharge(-1)
        mol.GetAtomWithIdx(oh).SetNumExplicitHs(0)
        mol.GetAtomWithIdx(oh).SetNoImplicit(True)
    out = mol.GetMol()
    Chem.SanitizeMol(out)
    return out, len(hits)


def build_ligand(pose_path, anion=False, resname="LIG"):
    """RDKit molecule with hydrogens at the docked position; returns (mol, n_acids_deprotonated, cip)."""
    heavy = pose_to_mol(pose_path)                       # heavy atoms, docked coordinates
    Chem.AssignStereochemistryFrom3D(heavy)
    centres = Chem.FindMolChiralCenters(heavy, useLegacyImplementation=False)
    cip = centres[0][1] if len(centres) == 1 else None
    heavy = Chem.RWMol(heavy)
    for a in heavy.GetAtoms():                           # RemoveHs leaves counts that block re-adding
        a.SetNumExplicitHs(0)
        a.SetNoImplicit(False)
    heavy = heavy.GetMol()
    Chem.SanitizeMol(heavy)
    n_acid = 0
    if anion:
        heavy, n_acid = deprotonate_acids(heavy)
    mol = Chem.AddHs(heavy, addCoords=True)
    for i, a in enumerate(mol.GetAtoms()):
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(resname)
        info.SetName(f"{a.GetSymbol():>2s}{i + 1:<2d}"[:4] if a.GetSymbol() != "H" else f"H{i + 1}"[:4])
        info.SetResidueNumber(1)
        info.SetChainId("L")
        info.SetIsHeteroAtom(True)
        a.SetMonomerInfo(info)
    return mol, n_acid, cip


def heavy_coords(mol):
    conf = mol.GetConformer()
    return [tuple(conf.GetAtomPosition(a.GetIdx())) for a in mol.GetAtoms() if a.GetAtomicNum() > 1]


def check_against_pose(mol, pose_path, cip):
    """Heavy-atom coordinates must be those of the pose and the stereocentre must be unchanged."""
    ref = heavy_coords(pose_to_mol(pose_path))
    got = heavy_coords(mol)
    if len(ref) != len(got) or max(max(abs(a - b) for a, b in zip(p, q)) for p, q in zip(ref, got)) > 1e-3:
        raise RuntimeError("heavy-atom coordinates changed during export")
    probe = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(probe)
    centres = Chem.FindMolChiralCenters(probe, useLegacyImplementation=False)
    if cip is not None and (len(centres) != 1 or centres[0][1] != cip):
        raise RuntimeError(f"stereocentre changed during export: {cip} -> {centres}")


def write_complex(site_pdb, ligand_mol, path):
    with open(site_pdb) as fh:
        site = [l for l in fh.read().splitlines() if l.startswith(("ATOM", "HETATM"))]
    block = [l for l in Chem.MolToPDBBlock(ligand_mol).splitlines() if l.startswith(("ATOM", "HETATM"))]
    with open(path, "w") as fh:
        fh.write("\n".join(site + ["TER"] + [l.replace("ATOM  ", "HETATM", 1) for l in block] + ["END"]) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pose", required=True, help="docked pose (.pdbqt)")
    p.add_argument("--site", required=True, help="1SA0_site.pdb written by dock.py")
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--anion", action="store_true", help="deprotonate carboxylic acids")
    p.add_argument("--resname", default="LIG")
    a = p.parse_args(argv)

    mol, n_acid, cip = build_ligand(a.pose, a.anion, a.resname)
    check_against_pose(mol, a.pose, cip)
    n_h = sum(1 for x in mol.GetAtoms() if x.GetAtomicNum() == 1)
    charge = Chem.GetFormalCharge(mol)
    w = Chem.SDWriter(a.out + ".sdf")
    w.write(mol)
    w.close()
    write_complex(a.site, mol, a.out + "_complex.pdb")
    print(f"ligand: {mol.GetNumAtoms()} atoms ({n_h} H), formal charge {charge:+d}, "
          f"{n_acid} carboxylic acid(s) deprotonated, CIP {cip}")
    print(f"wrote {a.out}.sdf and {a.out}_complex.pdb")


if __name__ == "__main__":
    main()
