\
#!/usr/bin/env python3
import argparse
from pathlib import Path

import requests
from Bio.PDB import PDBParser, PDBIO, Select

WATER_NAMES = {"HOH","WAT"}
COMMON_IONS = {"NA","CL","K","CA","MG","ZN","MN","FE","CU","SO4","PO4"}

class ProteinSelect(Select):
    def accept_atom(self, atom):
        return 1 if atom.get_parent().id[0] == " " else 0

class LigandSelect(Select):
    def __init__(self):
        self.best_res = None

    def pick_largest_ligand(self, structure):
        candidates = []
        for model in structure:
            for chain in model:
                for res in chain:
                    hetflag = res.id[0]
                    if hetflag == " ":
                        continue
                    resname = res.get_resname().strip()
                    if resname in WATER_NAMES or resname in COMMON_IONS:
                        continue
                    nat = len(list(res.get_atoms()))
                    candidates.append((nat, chain.id, res.id, resname))
        candidates.sort(reverse=True)
        return candidates[0] if candidates else None

    def accept_atom(self, atom):
        if self.best_res is None:
            return 0
        nat, chain_id, resid, resname = self.best_res
        res = atom.get_parent()
        chain = res.get_parent()
        return 1 if (chain.id == chain_id and res.id == resid) else 0

def run(cmd):
    import subprocess
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{p.stdout}")
    return p.stdout

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-id", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--obabel", default="obabel")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pdb_id = args.pdb_id.upper()
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_path = outdir/f"{pdb_id}.pdb"
    r = requests.get(pdb_url, timeout=60)
    r.raise_for_status()
    pdb_path.write_text(r.text)
    print("[OK] Saved:", pdb_path)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(pdb_path))

    io = PDBIO()
    io.set_structure(structure)

    receptor_pdb = outdir/f"{pdb_id}_protein.pdb"
    io.save(str(receptor_pdb), ProteinSelect())
    print("[OK] Receptor protein:", receptor_pdb)

    lig_sel = LigandSelect()
    lig_sel.best_res = lig_sel.pick_largest_ligand(structure)
    if lig_sel.best_res is None:
        raise RuntimeError("Could not find a suitable non-water, non-ion ligand to autobox.")
    print("[INFO] Picked ligand:", lig_sel.best_res)

    ref_lig_pdb = outdir/f"{pdb_id}_ref_lig.pdb"
    io.save(str(ref_lig_pdb), lig_sel)
    print("[OK] Ref ligand:", ref_lig_pdb)

    receptor_pdbqt = outdir/f"{pdb_id}_protein.pdbqt"
    ref_lig_sdf = outdir/f"{pdb_id}_ref_lig.sdf"

    # Convert using Open Babel (must be installed and on PATH)
    run([args.obabel, "-ipdb", str(receptor_pdb), "-opdbqt", "-O", str(receptor_pdbqt), "-xh"])
    run([args.obabel, "-ipdb", str(ref_lig_pdb), "-osdf", "-O", str(ref_lig_sdf), "-h"])

    print("[OK] Receptor PDBQT:", receptor_pdbqt)
    print("[OK] Ref ligand SDF:", ref_lig_sdf)

if __name__ == "__main__":
    main()
