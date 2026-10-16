"""Tests for the helper functions of docking/dock.py (no docking run needed).

Run:  python -m pytest tests/   (or python tests/test_docking_helpers.py)

Set SALMA_TEST_PDB to a local copy of 1SA0.pdb to also test the receptor-site construction.
"""
import math
import os
import sys
import tempfile

from rdkit import Chem

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "docking"))
import dock as D  # noqa: E402

CN2_LINES = open(os.path.join(HERE, "data", "CN2_1SA0_B700.pdb")).read().splitlines()
COL = D.BENCHMARKS["Colchicine"][0]


def coords(m):
    return [tuple(round(v, 4) for v in m.GetConformer().GetAtomPosition(i)) for i in range(m.GetNumAtoms())]


def test_all_ligand_definitions_are_valid():
    for n, s in D.SALMA.items():
        D.check_smiles(n, s)
    for n, (s, f) in D.BENCHMARKS.items():
        D.check_smiles(n, s, f)
    try:
        D.check_smiles("bad", "CCO", "C9H9")
    except ValueError:
        pass
    else:
        raise AssertionError("a wrong formula must be rejected")


def test_enantiomer_is_fixed_for_every_salma_compound():
    for name, smi in D.SALMA.items():
        mol = Chem.MolFromSmiles(smi)
        for want in ("S", "R"):
            m, k = D.embed_with_cip(mol, want)
            assert D.cip_of_single_centre(m) == want, (name, want)
            m2, k2 = D.embed_with_cip(mol, want)
            assert k == k2 and coords(m) == coords(m2), (name, want)    # deterministic


def test_embed_with_cip_needs_exactly_one_centre():
    for smi in ("CCO", D.BENCHMARKS["Podophyllotoxin"][0]):
        try:
            D.embed_with_cip(Chem.MolFromSmiles(smi), "S")
        except ValueError:
            pass
        else:
            raise AssertionError(smi)


def test_unconstrained_embedding_depends_on_rdkit_for_the_enantiomer():
    """Documents why the CIP label is enforced: over many seeds both enantiomers occur."""
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(D.SALMA["Salma-12"]))
    seen = set()
    for seed in range(40):
        m = Chem.Mol(mol)
        AllChem.EmbedMolecule(m, randomSeed=seed)
        seen.add(D.cip_of_single_centre(m))
    assert seen == {"R", "S"}


def test_crystal_cn2_geometry():
    cn2 = D.cn2_from_crystal(CN2_LINES)
    assert cn2.GetNumHeavyAtoms() == 30
    assert D.cip_of_single_centre(cn2) == "S"
    assert abs(D.biaryl_torsion(cn2) - 55.2) < 0.1


def test_colchicine_twist_is_controlled():
    mol = Chem.MolFromSmiles(COL)
    for sign in (1, -1):
        m, k, t = D.embed_with_twist(mol, sign)
        assert t * sign > 30 and abs(D.biaryl_torsion(m) - t) < 1e-6
        assert D.cip_of_single_centre(m) == "S"                          # C7 configuration unchanged
        m2, k2, t2 = D.embed_with_twist(mol, sign)
        assert (k, round(t, 3)) == (k2, round(t2, 3))
    # the two senses are different diastereomers in the pocket: torsions differ in sign, same magnitude
    assert abs(abs(D.embed_with_twist(mol, 1)[2]) - abs(D.embed_with_twist(mol, -1)[2])) < 5


def test_pose_from_colab_run_has_the_opposite_twist():
    """Colchicine_s1.pdbqt (first Colab run): twist -59 deg, C7 S; the crystal CN2 has +55 deg."""
    pose = D.pose_to_mol(os.path.join(HERE, "data", "Colchicine_s1_colab.pdbqt"))
    assert abs(D.biaryl_torsion(pose) - (-59.0)) < 0.2
    assert D.cip_of_single_centre(pose) == "S"


def test_centroid_of_cn2_is_the_published_box_centre():
    centre, n = D.ligand_centroid(CN2_LINES, *D.REF_LIGAND)
    assert n == 30
    assert tuple(round(c, 3) for c in centre) == D.PUBLISHED_CENTRE
    try:
        D.ligand_centroid(CN2_LINES, "XXX", "B", "1")
    except ValueError:
        pass
    else:
        raise AssertionError


def test_box_centre_check():
    D.check_box_centre(D.PUBLISHED_CENTRE)
    D.check_box_centre((117.2194, 90.1796, 6.2904))                      # rounding differences are fine
    for off in ((0.02, 0, 0), (0, -0.02, 0), (0, 0, 0.5)):
        shifted = tuple(a + b for a, b in zip(D.PUBLISHED_CENTRE, off))
        try:
            D.check_box_centre(shifted)
        except ValueError:
            pass
        else:
            raise AssertionError(off)


def test_box_distance_and_clearance():
    c = (10.0, 20.0, 30.0)
    assert D.dist_to_box((10, 20, 30), c, 22.0) == 0.0
    assert D.dist_to_box((21, 20, 30), c, 22.0) == 0.0                   # on the face
    assert abs(D.dist_to_box((24, 20, 30), c, 22.0) - 3.0) < 1e-12
    assert abs(D.dist_to_box((24, 25, 30), c, 22.0) - 3.0) < 1e-12       # still outside along x only
    assert abs(D.dist_to_box((24, 35, 30), c, 22.0) - 5.0) < 1e-12       # 3 in x, 4 in y
    line = lambda ch, resi, x, y, z: f"ATOM  {1:5d}  CA  ALA {ch}{resi:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
    pdb = [line("A", 5, 24, 20, 30), line("B", 7, 10, 20, 52), line("A", 9, 10, 20, 31)]
    d, where = D.deleted_residue_clearance(pdb, ["A:5", "B:7"], c, 22.0)
    assert where == "A:5" and abs(d - 3.0) < 1e-9                        # A:9 is not deleted, so it is ignored
    d, where = D.deleted_residue_clearance(pdb, ["A:9"], c, 22.0)
    assert d == 0.0                                                      # inside the box
    assert D.deleted_residue_clearance(pdb, [], c, 22.0) == (math.inf, None)
    assert D.deleted_residue_atoms(pdb, ["B:7"]) == [(10.0, 20.0, 52.0)]


def test_min_distance_and_pdbqt_xyz():
    assert D.min_distance([(0, 0, 0)], [(3, 4, 0), (1, 0, 0)]) == 1.0
    assert D.min_distance([], [(1, 0, 0)]) == math.inf
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "x.pdbqt")
        open(p, "w").write("MODEL 1\nATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C \n"
                           "ENDMDL\nMODEL 2\nATOM      1  C   UNL     1       9.000   9.000   9.000  1.00  0.00    +0.000 C \nENDMDL\n")
        assert D.pdbqt_xyz(p) == [(1.0, 2.0, 3.0)]                       # first model only


def test_removed_residue_parser_on_the_real_meeko_message():
    log = ("noise\n- Template matching failed for: ['A:35', 'A:47', 'B:437'] Ignored due to allow_bad_res.\nmore\n")
    assert D.deleted_residues(log) == ["A:35", "A:47", "B:437"]
    assert D.deleted_residues("nothing removed") == []
    try:
        D.deleted_residues("No template matched for residue_key='A:35'")
    except RuntimeError:
        pass
    else:
        raise AssertionError("an unparseable list must not pass silently")


def test_names_and_statistics():
    assert D.tag_of("Combretastatin A-4") == "Combretastatin_A-4" and D.tag_of("CN2 (re-dock)") == "CN2_re-dock_"
    mu, sd = D.mean_sd([1.0, 2.0, 3.0])
    assert mu == 2.0 and abs(sd - 1.0) < 1e-12                           # sample standard deviation
    assert D.mean_sd([5.0]) == (5.0, 0.0)


def test_unknown_ligand_name_is_rejected_before_any_work():
    with tempfile.TemporaryDirectory() as d:
        try:
            D.main(["--workdir", d, "--only", "Salma-99"])
        except SystemExit as e:
            assert e.code == 2
        else:
            raise AssertionError


def test_site_construction_on_the_full_structure():
    path = os.environ.get("SALMA_TEST_PDB")
    if not path or not os.path.exists(path):
        print("  (skipped: set SALMA_TEST_PDB to a local 1SA0.pdb)")
        return
    lines = D.read_lines(path)
    site = D.build_site(lines)
    assert (sum(l.startswith("ATOM") for l in site), sum(l.startswith("HETATM") for l in site)) == (6527, 62)
    centre, n = D.ligand_centroid(lines, *D.REF_LIGAND)
    assert n == 30 and all(abs(a - b) < 0.001 for a, b in zip(centre, D.PUBLISHED_CENTRE))
    try:
        D.build_site(lines[:-2000])
    except ValueError:
        pass
    else:
        raise AssertionError("a truncated structure must be rejected")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print("PASS", t.__name__)
    print(f"{len(tests)} tests passed")
