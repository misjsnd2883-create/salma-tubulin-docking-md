#!/usr/bin/env python3
"""Docking of the Salma series and reference ligands into the colchicine site of tubulin.

Protocol (manuscript Methods 2.3):
  * receptor: PDB 1SA0, chains A+B; GTP (A600), Mg (A501), GDP (B602) and Mg (B502) are kept
  * search box: 22 x 22 x 22 A, centred on the centroid of the crystallographic ligand CN2
    (chain B, residue 700), computed here from the deposited coordinates
  * AutoDock Vina 1.2.7, exhaustiveness 32, five seeds (1-5); Meeko 0.7.1 (`-a`)
  * ligands are embedded with RDKit (seed 0xf00d) and optimised with MMFF

Choices that would otherwise be left to chance are fixed explicitly:
  * the Salma compounds have one stereocentre that is undefined in the SMILES; the embedding
    is repeated (seed 0xf00d + k) until the centre has the requested CIP label (default S)
  * colchicine is chiral about its biaryl axis; the embedding is repeated until the biaryl
    twist has the same sense as the crystallographic CN2 (`--colchicine-twist crystal`,
    default); `--colchicine-twist free` takes the first embedding, as an unconstrained run would

Validation: CN2 is re-docked from its crystal geometry and the heavy-atom RMSD to the crystal
pose is reported (`--no-redock` to skip). Output: <workdir>/results.json plus one best pose
(.pdbqt) per ligand and seed.

Usage:
    python dock.py --workdir work
    python dock.py --workdir work --only Salma-12 Colchicine --seeds 1 2
"""
import argparse
import glob
import hashlib
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
import time
import urllib.request
from datetime import datetime, timezone
from importlib import metadata

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import rdMolTransforms as T
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# ----------------------------------------------------------------------------- settings
PDB_ID = "1SA0"
KEEP_CHAINS = ("A", "B")
KEEP_HET = {("GTP", "A", "600"), ("MG", "A", "501"), ("GDP", "B", "602"), ("MG", "B", "502")}
REF_LIGAND = ("CN2", "B", "700")
EXPECTED_SITE_ATOMS = (6527, 62)          # ATOM, HETATM lines of chains A+B plus the kept cofactors
PUBLISHED_CENTRE = (117.219, 90.180, 6.290)
BOX = 22.0
EXHAUSTIVENESS = 32
SEEDS = [1, 2, 3, 4, 5]
EMBED_SEED = 0xF00D
VINA_CUTOFF = 8.0                          # A; Vina's interaction cutoff

SALMA = {
    "Salma-1": "COc1cc(F)c(CC(C)C(=O)O)cc1O",
    "Salma-3": "COc1cc(F)c(CC(C(C)C)C(=O)O)cc1F",
    "Salma-4": "COc1cc(F)c(CC(C(C)C)C(=O)NC)cc1F",
    "Salma-5": "Oc1cc(F)c(CC(C(C)C)C(=O)O)cc1F",
    "Salma-6": "COc1cc(F)c(CC(C2CC2)C(=O)O)cc1F",
    "Salma-12": "COc1c(C(=O)N)c(F)c(CC(C2CCCCC2)C(=O)O)cc1F",
    "Salma-20": "COc1c(F)cc(CC(C(=O)O)C(C)(C)C)c(F)c1C#N",
}
BENCHMARKS = {  # name: (SMILES, expected molecular formula)
    "Colchicine": ("CC(=O)N[C@H]1CCc2cc(OC)c(OC)c(OC)c2-c2ccc(OC)c(=O)cc21", "C22H25NO6"),
    "Combretastatin A-4": (r"COc1ccc(/C=C\c2cc(OC)c(OC)c(OC)c2)cc1O", "C18H20O5"),
    "Podophyllotoxin": ("COc1cc(cc(OC)c1OC)[C@@H]1[C@@H]2[C@H](COC2=O)[C@H](O)c2cc3OCOc3cc12", "C22H22O8"),
    "Nocodazole": ("COC(=O)Nc1nc2ccc(C(=O)c3cccs3)cc2[nH]1", "C14H11N3O3S"),
}
# CN2 = DAMA-colchicine; template used to assign bond orders to the crystal coordinates
CN2_TEMPLATE = "COc1cc2c(c(OC)c1OC)-c1ccc(OC)c(=O)cc1C(NC(=O)CS)CC2"
BIARYL = Chem.MolFromSmarts("[CH2]-c:c-c:c-[CH1]")   # torsion across the biaryl axis of ring B


def tag_of(name):
    """File-name-safe version of a ligand name."""
    return re.sub(r"[^A-Za-z0-9_-]+", "_", name)


# ------------------------------------------------------------------------ PDB handling
def fetch_pdb(dest, pdb_id=PDB_ID):
    urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.pdb", dest)


def read_lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


def build_site(lines):
    """Chains A+B plus the retained cofactors; the atom counts are checked against 1SA0."""
    out = []
    for line in lines:
        rec = line[:6].strip()
        if rec == "ATOM" and line[21] in KEEP_CHAINS:
            out.append(line)
        elif rec == "HETATM" and (line[17:20].strip(), line[21], line[22:26].strip()) in KEEP_HET:
            out.append(line)
    counts = (sum(1 for l in out if l.startswith("ATOM")), sum(1 for l in out if l.startswith("HETATM")))
    if counts != EXPECTED_SITE_ATOMS:
        raise ValueError(f"unexpected atom counts {counts}, expected {EXPECTED_SITE_ATOMS} for {PDB_ID}")
    return out


def residue_lines(lines, resname, chain, resid):
    return [l for l in lines if l.startswith("HETATM") and l[17:20].strip() == resname
            and l[21] == chain and l[22:26].strip() == resid]


def ligand_centroid(lines, resname, chain, resid):
    """Centroid of the deposited atoms of one hetero residue; returns ((x, y, z), n_atoms)."""
    xyz = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in residue_lines(lines, resname, chain, resid)]
    if not xyz:
        raise ValueError(f"{resname} {chain}{resid} not found")
    n = len(xyz)
    return tuple(sum(p[i] for p in xyz) / n for i in range(3)), n


def check_box_centre(centre, published=PUBLISHED_CENTRE, tol=0.01):
    """The centroid computed from the deposited coordinates must agree with the centre used in the paper."""
    if any(abs(c - e) > tol for c, e in zip(centre, published)):
        raise ValueError(f"centroid of {REF_LIGAND[0]} {tuple(centre)} differs from the published box centre {tuple(published)}")


# ------------------------------------------------------------------------------ receptor
_PREP_PROGRAM = r"""
import os, sys
ccd = os.environ.get("SALMA_CCD_DIR")
if ccd:   # offline use: read the GDP/GTP chemical-component files from a local directory
    import meeko.chemtempgen as ct
    ct.fetch_from_pdb = lambda resname, *a, **k: os.path.join(ccd, resname + ".cif")
from meeko.cli.mk_prepare_receptor import main
sys.argv = ["mk_prepare_receptor.py"] + sys.argv[1:]
main()
"""


def prepare_receptor(site_pdb, out_pdbqt, centre, box, log_path, ccd_dir=None):
    """Run Meeko's receptor preparation (`mk_prepare_receptor.py ... -a`). Needs internet access
    to download the GDP/GTP templates from RCSB unless `ccd_dir` holds GDP.cif and GTP.cif."""
    env = dict(os.environ)
    if ccd_dir:
        env["SALMA_CCD_DIR"] = os.path.abspath(ccd_dir)
    cmd = [sys.executable, "-c", _PREP_PROGRAM, "--read_pdb", site_pdb, "-p", out_pdbqt, "-a",
           "--box_center", *[str(c) for c in centre], "--box_size", *[str(box)] * 3]
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    text = res.stdout + res.stderr
    with open(log_path, "w") as fh:
        fh.write(text)
    if res.returncode != 0 or not os.path.exists(out_pdbqt):
        raise RuntimeError(f"receptor preparation failed (see {log_path}):\n{text[-1500:]}")
    return text


def deleted_residues(log_text):
    """Residues that Meeko could not match to a template and therefore removed (e.g. 'B:326')."""
    m = re.search(r"Template matching failed for: \[(.*?)\]", log_text, re.S)
    if not m:
        if "No template matched" in log_text:
            raise RuntimeError("Meeko reported unmatched residues but the list could not be parsed")
        return []
    return re.findall(r"'([A-Z]:\d+)'", m.group(1))


def dist_to_box(p, centre, box):
    """Euclidean distance (A) from point p to the search cube; 0 if p is inside."""
    half = box / 2.0
    return math.sqrt(sum(max(abs(p[i] - centre[i]) - half, 0.0) ** 2 for i in range(3)))


def deleted_residue_clearance(pdb_lines, deleted, centre, box):
    """Smallest distance from the search box to any atom of a deleted residue.

    Returns (distance in A, residue key); (inf, None) if nothing was deleted."""
    wanted = set(deleted)
    best, where = math.inf, None
    for l in pdb_lines:
        if l[:4] != "ATOM":
            continue
        key = f"{l[21]}:{l[22:26].strip()}"
        if key in wanted:
            d = dist_to_box((float(l[30:38]), float(l[38:46]), float(l[46:54])), centre, box)
            if d < best:
                best, where = d, key
    return best, where


def deleted_residue_atoms(pdb_lines, deleted):
    """Coordinates of the deposited atoms of the residues that Meeko removed."""
    wanted = set(deleted)
    return [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in pdb_lines
            if l[:4] == "ATOM" and f"{l[21]}:{l[22:26].strip()}" in wanted]


def pdbqt_xyz(path):
    """Atom coordinates of the first model of a PDBQT file."""
    out = []
    with open(path) as fh:
        for l in fh:
            if l.startswith("ENDMDL"):
                break
            if l.startswith(("ATOM", "HETATM")):
                out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def min_distance(points_a, points_b):
    """Smallest distance (A) between two sets of points; inf if either is empty."""
    return min((math.dist(p, q) for p in points_a for q in points_b), default=math.inf)


# --------------------------------------------------------------------------------- maps
def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def build_maps(receptor, prefix, centre, box, stamp_path):
    """Compute the Vina maps once; they are reused only if receptor, centre and box are unchanged."""
    from vina import Vina
    stamp = json.dumps({"receptor": _sha(receptor), "centre": list(centre), "box": box})
    if glob.glob(prefix + ".*.map") and os.path.exists(stamp_path) and open(stamp_path).read() == stamp:
        return "reused"
    for f in glob.glob(prefix + ".*"):
        if f.endswith((".map", ".fld", ".xyz")):
            os.remove(f)
    v = Vina(sf_name="vina", seed=1, verbosity=0)
    v.set_receptor(receptor)
    v.compute_vina_maps(center=list(centre), box_size=[box] * 3, force_even_voxels=True)
    v.write_maps(prefix)
    with open(stamp_path, "w") as fh:
        fh.write(stamp)
    return "built"


# ----------------------------------------------------------------------------- ligands
def check_smiles(name, smiles, formula=None):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"{name}: SMILES does not parse")
    if not any(len(r) == 6 and all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r) for r in m.GetRingInfo().AtomRings()):
        raise ValueError(f"{name}: no six-membered aromatic ring")
    if formula and CalcMolFormula(m) != formula:
        raise ValueError(f"{name}: formula {CalcMolFormula(m)} != expected {formula}")
    return m


def _mmff(m):
    """MMFF optimisation; returns True if parameters exist (convergence is not required)."""
    return AllChem.MMFFOptimizeMolecule(m, maxIters=500) in (0, 1)


def embed_plain(mol, embed_seed=EMBED_SEED):
    """First embedding at the given seed, MMFF-optimised (the unconstrained behaviour)."""
    m = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(m, randomSeed=embed_seed) != 0:
        raise RuntimeError("embedding failed")
    if not _mmff(m):
        raise RuntimeError("MMFF parameters missing")
    return m, 0


def cip_of_single_centre(m):
    Chem.AssignStereochemistryFrom3D(m)
    c = Chem.FindMolChiralCenters(Chem.RemoveHs(m), useLegacyImplementation=False)
    return c[0][1] if len(c) == 1 else None


def embed_with_cip(mol, want="S", embed_seed=EMBED_SEED, max_tries=200):
    """Embed until the single stereocentre has the requested CIP label.

    RDKit chooses the enantiomer at random when the SMILES leaves the centre undefined, so a
    fixed seed alone does not fix the enantiomer across RDKit versions. Returns (mol, k) where
    k is the number of rejected embeddings."""
    m0 = Chem.AddHs(mol)
    if len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)) != 1:
        raise ValueError("expected exactly one stereocentre")
    for k in range(max_tries):
        m = Chem.Mol(m0)
        if AllChem.EmbedMolecule(m, randomSeed=embed_seed + k) != 0:
            continue
        if cip_of_single_centre(m) == want and _mmff(m):
            return m, k
    raise RuntimeError(f"no {want} embedding found in {max_tries} tries")


def biaryl_torsion(m):
    """Torsion (deg) across the biaryl axis of colchicine-type ligands."""
    match = m.GetSubstructMatches(BIARYL)
    if len(match) != 1:
        raise ValueError(f"expected one biaryl match, found {len(match)}")
    a = match[0]
    return T.GetDihedralDeg(m.GetConformer(), a[1], a[2], a[3], a[4])


def embed_with_twist(mol, want_sign, embed_seed=EMBED_SEED, max_tries=200, min_abs=30.0):
    """Embed until the biaryl torsion (after MMFF) has the requested sign.

    Returns (mol, k, torsion)."""
    m0 = Chem.AddHs(mol)
    for k in range(max_tries):
        m = Chem.Mol(m0)
        if AllChem.EmbedMolecule(m, randomSeed=embed_seed + k) != 0 or not _mmff(m):
            continue
        t = biaryl_torsion(m)
        if t * want_sign > min_abs:
            return m, k, t
    raise RuntimeError("no embedding with the requested biaryl twist found")


def to_pdbqt(mol_h):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    s, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mol_h)[0])
    if not ok:
        raise RuntimeError(err)
    return s


def pose_to_mol(path):
    """Docked pose (first model of a PDBQT file) as an RDKit molecule without hydrogens."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    mol = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(path, skip_typing=True))[0]
    return Chem.RemoveHs(mol)


def cn2_from_crystal(lines):
    """CN2 with bond orders assigned to the deposited coordinates (heavy atoms only)."""
    blk = "\n".join(residue_lines(lines, *REF_LIGAND))
    raw = Chem.MolFromPDBBlock(blk, removeHs=False, proximityBonding=True, sanitize=False)
    mol = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(CN2_TEMPLATE), raw)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    return mol


# ----------------------------------------------------------------------------- docking
def run_dock(ligand_pdbqt, seed, pose_path, maps_prefix, exhaustiveness, cpu=0):
    from vina import Vina
    v = Vina(sf_name="vina", seed=seed, verbosity=0, cpu=cpu)
    v.load_maps(maps_prefix)
    v.set_ligand_from_string(ligand_pdbqt)
    v.dock(exhaustiveness=exhaustiveness, n_poses=5)
    v.write_poses(pose_path, n_poses=1, overwrite=True)
    return float(v.energies(n_poses=1)[0][0])


def mean_sd(values):
    return statistics.mean(values), (statistics.stdev(values) if len(values) > 1 else 0.0)


def versions():
    out = {"python": platform.python_version(), "rdkit": rdBase.rdkitVersion}
    for pkg in ("meeko", "vina", "gemmi", "numpy"):
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = None
    return out


# ------------------------------------------------------------------------------- driver
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--workdir", default="docking_work")
    p.add_argument("--pdb", help="local copy of 1SA0.pdb (default: download from RCSB)")
    p.add_argument("--ccd-dir", help="directory with GDP.cif and GTP.cif (offline receptor preparation)")
    p.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    p.add_argument("--exhaustiveness", type=int, default=EXHAUSTIVENESS)
    p.add_argument("--only", nargs="+", help="dock only these ligands (names as printed in the table)")
    p.add_argument("--enantiomer", choices=("S", "R"), default="S", help="CIP label for the Salma compounds")
    p.add_argument("--colchicine-twist", choices=("crystal", "free"), default="crystal")
    p.add_argument("--no-redock", action="store_true", help="skip the CN2 re-docking validation")
    p.add_argument("--cpu", type=int, default=0, help="CPU cores for Vina (0 = all)")
    a = p.parse_args(argv)

    os.makedirs(a.workdir, exist_ok=True)
    W = lambda name: os.path.join(a.workdir, name)

    # ---- ligand definitions
    mols = {n: check_smiles(n, s) for n, s in SALMA.items()}
    mols.update({n: check_smiles(n, s, f) for n, (s, f) in BENCHMARKS.items()})
    names = list(mols)
    if a.only:
        by_key = {k.lower(): n for n in names for k in (n, tag_of(n))}
        unknown = [x for x in a.only if x.lower() not in by_key]
        if unknown:
            p.error(f"unknown ligand(s): {unknown}; choose from {names}")
        keep = {by_key[x.lower()] for x in a.only}
        names = [n for n in names if n in keep]
    print(f"{len(mols)} ligand definitions validated; docking {len(names)}"
          f"{' + CN2 re-dock' if not a.no_redock else ''}", flush=True)

    # ---- receptor
    pdb_path = a.pdb or W(f"{PDB_ID}.pdb")
    if not os.path.exists(pdb_path):
        fetch_pdb(pdb_path)
    lines = read_lines(pdb_path)
    site = build_site(lines)
    open(W("1SA0_site.pdb"), "w").writelines(l + "\n" for l in site + ["END"])
    centre, n_ref = ligand_centroid(lines, *REF_LIGAND)
    centre = tuple(round(c, 3) for c in centre)
    check_box_centre(centre)
    print(f"box centre = centroid of {n_ref} {REF_LIGAND[0]} atoms: {centre}  (published {PUBLISHED_CENTRE})", flush=True)

    receptor = W("receptor.pdbqt")
    log_path = W("prep_log.txt")
    if not os.path.exists(receptor) or not os.path.exists(log_path):
        prepare_receptor(W("1SA0_site.pdb"), receptor, centre, BOX, log_path, a.ccd_dir)
    log_text = open(log_path).read()
    deleted = deleted_residues(log_text)
    clearance, where = deleted_residue_clearance(lines, deleted, centre, BOX)
    removed_atoms = deleted_residue_atoms(lines, deleted)
    ref_xyz = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in residue_lines(lines, *REF_LIGAND)]
    print(f"residues removed by Meeko: {len(deleted)}; nearest to the box: {where} ({clearance:.1f} A from its edge); "
          f"nearest to the crystal {REF_LIGAND[0]}: {min_distance(removed_atoms, ref_xyz):.1f} A", flush=True)
    if clearance <= 0:
        raise RuntimeError(f"a removed residue ({where}) has atoms inside the search box")
    rec_txt = open(receptor).read()
    for cof in ("GTP", "GDP", " MG"):
        if cof not in rec_txt:
            raise RuntimeError(f"cofactor {cof.strip()} missing from the prepared receptor")
    print("maps:", build_maps(receptor, W("grid"), centre, BOX, W("maps.stamp")), flush=True)

    # ---- reference twist from the crystal ligand
    cn2 = cn2_from_crystal(lines)
    crystal_twist = biaryl_torsion(cn2)
    print(f"crystal CN2: biaryl torsion {crystal_twist:+.1f} deg", flush=True)

    # ---- results file (resumable)
    meta = {"pdb": PDB_ID, "centre": list(centre), "box": BOX, "exhaustiveness": a.exhaustiveness,
            "embed_seed": hex(EMBED_SEED), "enantiomer": a.enantiomer, "colchicine_twist": a.colchicine_twist,
            "receptor_sha256": _sha(receptor)}
    res_path = W("results.json")
    data = {"meta": {**meta, "versions": versions(), "created": datetime.now(timezone.utc).isoformat()},
            "scores": {}, "ligands": {}, "redock": {}, "removed_residue_distance": {}}
    if os.path.exists(res_path):
        old = json.load(open(res_path))
        if {k: old["meta"].get(k) for k in meta} != meta:
            raise RuntimeError(f"{res_path} was produced with different settings; use another --workdir")
        data = old
        data.setdefault("removed_residue_distance", {})
    save = lambda: json.dump(data, open(res_path, "w"), indent=1)

    def dock_set(name, pdbqt, info):
        tag = tag_of(name)
        data["ligands"][name] = info
        scores = data["scores"].setdefault(name, {})
        for seed in a.seeds:
            if str(seed) in scores:
                continue
            t0 = time.time()
            pose = W(f"{tag}_s{seed}.pdbqt")
            scores[str(seed)] = run_dock(pdbqt, seed, pose, W("grid"), a.exhaustiveness, a.cpu)
            near = min_distance(removed_atoms, pdbqt_xyz(pose))
            data["removed_residue_distance"].setdefault(name, {})[str(seed)] = round(near, 2)
            save()
            flag = f"  <-- within Vina's {VINA_CUTOFF:g} A cutoff of a removed residue" if near < VINA_CUTOFF else ""
            print(f"  {name:20s} seed {seed}: {scores[str(seed)]:7.3f}  ({time.time() - t0:.0f} s){flag}", flush=True)
        return tag

    # ---- CN2 re-docking validation
    if not a.no_redock:
        print("CN2 re-docking (crystal geometry as input)", flush=True)
        cn2_h = Chem.AddHs(cn2, addCoords=True)
        tag = dock_set("CN2 (re-dock)", to_pdbqt(cn2_h), {"input": "crystal geometry"})
        rms = data["redock"].setdefault("rmsd", {})
        for seed in a.seeds:
            if str(seed) not in rms:
                rms[str(seed)] = rdMolAlign.CalcRMS(pose_to_mol(W(f"{tag}_s{seed}.pdbqt")), Chem.RemoveHs(cn2))
        save()

    # ---- ligands
    for name in names:
        print(name, flush=True)
        mol = mols[name]
        info = {"heavy_atoms": mol.GetNumHeavyAtoms()}
        if name in SALMA:
            m, k = embed_with_cip(mol, a.enantiomer)
            info.update(cip=a.enantiomer, embed_rejections=k)
        elif name == "Colchicine" and a.colchicine_twist == "crystal":
            m, k, t = embed_with_twist(mol, 1 if crystal_twist > 0 else -1)
            info.update(embed_rejections=k, embedded_biaryl_torsion=round(t, 1))
        else:
            m, k = embed_plain(mol)
            info.update(embed_rejections=k)
        tag = dock_set(name, to_pdbqt(m), info)
        if name in SALMA:  # confirm the enantiomer survived docking
            c = cip_of_single_centre(pose_to_mol(W(f"{tag}_s{a.seeds[0]}.pdbqt")))
            if c != a.enantiomer:
                print(f"  WARNING: docked pose of {name} has CIP {c}, expected {a.enantiomer}", flush=True)

    # ---- summary
    print("\nligand               heavy   mean +- SD (kcal/mol)   LE (kcal/mol/heavy atom)   n   closest removed residue (A)")
    for name, sc in data["scores"].items():
        vals = list(sc.values())
        mu, sd = mean_sd(vals)
        heavy = data["ligands"].get(name, {}).get("heavy_atoms") or (n_ref if name.startswith("CN2") else None)
        le = f"{-mu / heavy:.3f}" if heavy else "-"
        near = min(data["removed_residue_distance"].get(name, {"": math.nan}).values())
        print(f"{name:20s} {heavy if heavy else '-':>5}   {mu:7.3f} +- {sd:5.3f}            {le:>8}              {len(vals)}   {near:.1f}")
    if data["redock"].get("rmsd"):
        mu, sd = mean_sd(list(data["redock"]["rmsd"].values()))
        print(f"\nCN2 re-dock heavy-atom RMSD to the crystal pose: {mu:.2f} +- {sd:.2f} A")
    print(f"\nresults: {res_path}")


if __name__ == "__main__":
    main()
