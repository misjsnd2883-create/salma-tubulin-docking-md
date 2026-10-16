"""Build GDP.cif and GTP.cif from the coordinates in 1SA0, for offline runs of docking/dock.py.

Meeko normally downloads the chemical-component definitions of GDP and GTP from RCSB. Where the
network blocks that, `dock.py --ccd-dir <dir>` can read local files instead. These are minimal
CIF files (atom names, elements, bond orders) assembled from the 1SA0 coordinates and the
neutral chemical structures of GDP and GTP; they are meant for testing and offline use, not as
replacements for the official definitions.

Usage:
    python make_ccd_from_crystal.py 1SA0.pdb ccd_dir
"""
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

SPEC = {  # resname: (chain, resid, neutral structure, stereo ignored)
    "GDP": ("B", "602", "Nc1nc2c(ncn2C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"),
    "GTP": ("A", "600", "Nc1nc2c(ncn2C2OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C2O)c(=O)[nH]1"),
}


def make_ccd(pdb_lines, outdir):
    os.makedirs(outdir, exist_ok=True)
    for res, (chain, resid, smiles) in SPEC.items():
        block = "\n".join(l for l in pdb_lines if l.startswith("HETATM") and l[17:20].strip() == res
                          and l[21] == chain and l[22:26].strip() == resid)
        raw = Chem.MolFromPDBBlock(block, removeHs=False, proximityBonding=True, sanitize=False)
        template = Chem.MolFromSmiles(smiles)
        if raw.GetNumAtoms() != template.GetNumAtoms():
            raise ValueError(f"{res}: {raw.GetNumAtoms()} atoms in the crystal, {template.GetNumAtoms()} in the template")
        mol = AllChem.AssignBondOrdersFromTemplate(template, raw)
        Chem.SanitizeMol(mol)
        names = [a.GetPDBResidueInfo().GetName().strip() for a in mol.GetAtoms()]
        if len(set(names)) != len(names):
            raise ValueError(f"{res}: duplicate atom names")
        mh = Chem.AddHs(mol)
        names += [f"HX{i + 1}" for i in range(mh.GetNumAtoms() - mol.GetNumAtoms())]
        q = lambda s: f'"{s}"' if "'" in s else s
        out = [f"data_{res}", "#", "loop_", "_chem_comp_atom.comp_id", "_chem_comp_atom.atom_id",
               "_chem_comp_atom.type_symbol", "_chem_comp_atom.charge", "_chem_comp_atom.pdbx_leaving_atom_flag"]
        out += [f"{res} {q(names[a.GetIdx()])} {a.GetSymbol().upper()} 0 N" for a in mh.GetAtoms()]
        out += ["#", "loop_", "_chem_comp_bond.comp_id", "_chem_comp_bond.atom_id_1",
                "_chem_comp_bond.atom_id_2", "_chem_comp_bond.value_order"]
        kek = Chem.Mol(mh)
        Chem.Kekulize(kek, clearAromaticFlags=True)
        order = {Chem.BondType.SINGLE: "SING", Chem.BondType.DOUBLE: "DOUB"}
        out += [f"{res} {q(names[b.GetBeginAtomIdx()])} {q(names[b.GetEndAtomIdx()])} {order[b.GetBondType()]}"
                for b in kek.GetBonds()]
        with open(os.path.join(outdir, f"{res}.cif"), "w") as fh:
            fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    with open(sys.argv[1]) as fh:
        make_ccd(fh.read().splitlines(), sys.argv[2])
    print("wrote", os.path.join(sys.argv[2], "GDP.cif"), "and GTP.cif")
