"""export_pose.py on a Salma-12 pose written by dock.py (Vina 1.2.7, seed 2)."""
import os
import sys
import tempfile

from rdkit import Chem

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "docking"))
import export_pose as E  # noqa: E402

POSE = os.path.join(HERE, "data", "Salma-12_s2.pdbqt")
SITE = os.path.join(HERE, "data", "CN2_1SA0_B700.pdb")      # any PDB site file works for the complex writer


def test_anion_export():
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "S12")
        E.main(["--pose", POSE, "--site", SITE, "--out", out, "--anion"])
        mol = next(Chem.SDMolSupplier(out + ".sdf", removeHs=False))
        assert mol.GetNumAtoms() == 44
        assert sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1) == 20
        assert Chem.GetFormalCharge(mol) == -1
        Chem.AssignStereochemistryFrom3D(mol)
        assert [c[1] for c in Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)] == ["S"]
        lines = open(out + "_complex.pdb").read().splitlines()
        assert sum(1 for l in lines if l[17:20] == "LIG") == 44 and lines[-1] == "END"


def test_neutral_export_has_acid_hydrogen():
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "S12n")
        E.main(["--pose", POSE, "--site", SITE, "--out", out])
        mol = next(Chem.SDMolSupplier(out + ".sdf", removeHs=False))
        assert Chem.GetFormalCharge(mol) == 0 and mol.GetNumAtoms() == 45
