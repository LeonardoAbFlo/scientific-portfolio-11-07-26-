\
#!/usr/bin/env python3
import argparse, json, math, os, subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors

PHYS_DESC_FUNCS = {
    "MolWt": Descriptors.MolWt,
    "LogP": Crippen.MolLogP,
    "TPSA": rdMolDescriptors.CalcTPSA,
    "HBD": Lipinski.NumHDonors,
    "HBA": Lipinski.NumHAcceptors,
    "RotB": Lipinski.NumRotatableBonds,
    "RingCount": Lipinski.RingCount,
    "AromaticRings": Lipinski.NumAromaticRings,
    "AliphaticRings": Lipinski.NumAliphaticRings,
    "FractionCSP3": Lipinski.FractionCSP3,
    "HeavyAtomCount": Lipinski.HeavyAtomCount,
    "NHOHCount": Lipinski.NHOHCount,
    "NOCount": Lipinski.NOCount,
    "MolMR": Crippen.MolMR,
}

def canonicalize_smiles(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return Chem.MolToSmiles(m, isomericSmiles=True)

def compute_physchem(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return {k: float(fn(m)) for k, fn in PHYS_DESC_FUNCS.items()}

def smiles_to_3d_sdf(smiles, out_sdf):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    ok = AllChem.EmbedMolecule(m, params)
    if ok != 0:
        return False
    AllChem.UFFOptimizeMolecule(m, maxIters=200)
    w = Chem.SDWriter(str(out_sdf))
    w.write(m)
    w.close()
    return True

def run_gnina(gnina_bin, receptor_pdbqt, ligand_sdf, ref_lig_sdf, out_sdf, log_txt,
             exhaustiveness=8, num_modes=9, seed=0):
    cmd = [
        str(gnina_bin),
        "-r", str(receptor_pdbqt),
        "-l", str(ligand_sdf),
        "--autobox_ligand", str(ref_lig_sdf),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "--out", str(out_sdf),
        "--log", str(log_txt),
    ]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return p.returncode, p.stdout

def parse_gnina_sdf_best(out_sdf):
    suppl = Chem.SDMolSupplier(str(out_sdf), removeHs=False)
    best = None
    for mol in suppl:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()

        def get_prop(name):
            for key in props.keys():
                if key.lower() == name.lower():
                    return props[key]
            return None

        ma = get_prop("minimizedAffinity")
        cs = get_prop("CNNscore")
        ca = get_prop("CNNaffinity")

        try:
            ma = float(ma) if ma is not None else None
            cs = float(cs) if cs is not None else None
            ca = float(ca) if ca is not None else None
        except Exception:
            continue

        if ma is None:
            continue
        if (best is None) or (ma < best["minimizedAffinity"]):
            best = {"minimizedAffinity": ma, "CNNscore": cs, "CNNaffinity": ca}
    return best

def dock_one(row, receptor_pdbqt, ref_lig_sdf, gnina_bin, lig3d_dir, out_dir,
             exhaustiveness, num_modes, seed):
    mid = row["molecule_chembl_id"]
    smi = row["smiles_can"]

    lig_sdf = lig3d_dir / f"{mid}.sdf"
    if not lig_sdf.exists():
        ok = smiles_to_3d_sdf(smi, lig_sdf)
        if not ok:
            return None

    out_sdf = out_dir / f"{mid}_docked.sdf"
    log_txt = out_dir / f"{mid}.log"

    # resume: if already docked, reuse
    if out_sdf.exists():
        best = parse_gnina_sdf_best(out_sdf)
        if best is None:
            return None
    else:
        rc, stdout = run_gnina(
            gnina_bin=gnina_bin,
            receptor_pdbqt=receptor_pdbqt,
            ligand_sdf=lig_sdf,
            ref_lig_sdf=ref_lig_sdf,
            out_sdf=out_sdf,
            log_txt=log_txt,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            seed=seed,
        )
        if rc != 0 or (not out_sdf.exists()):
            log_txt.write_text(stdout)
            return None
        best = parse_gnina_sdf_best(out_sdf)
        if best is None:
            return None

    desc = compute_physchem(smi)
    if desc is None:
        return None

    return {
        "molecule_chembl_id": mid,
        "smiles": smi,
        "standard_value_nM_median": float(row["standard_value_nM_median"]),
        "n_activity_points": int(row["n_activity_points"]),
        "pchembl_median": float(row["pchembl_median"]) if pd.notnull(row["pchembl_median"]) else np.nan,
        **desc,
        **best
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--agg-csv", required=True, help="outputs/chembl/gba1_agg_medians.csv")
    ap.add_argument("--receptor-pdbqt", required=True)
    ap.add_argument("--ref-lig-sdf", required=True)
    ap.add_argument("--gnina", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-ligands", type=int, default=163)
    ap.add_argument("--exhaustiveness", type=int, default=64)
    ap.add_argument("--num-modes", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=1)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    dataset_dir = outdir/"datasets"
    dock_dir = outdir/"docking"
    lig3d_dir = dock_dir/"ligands_3d"
    gnina_out = dock_dir/"gnina_out"
    for p in [dataset_dir, lig3d_dir, gnina_out]:
        p.mkdir(parents=True, exist_ok=True)

    df_base = pd.read_csv(args.agg_csv)
    df_base["smiles_can"] = df_base["canonical_smiles"].apply(canonicalize_smiles)
    df_base = df_base.dropna(subset=["smiles_can"]).drop_duplicates("smiles_can").reset_index(drop=True)
    df_work = df_base.head(args.n_ligands).copy()
    print(f"[INFO] Docking N={len(df_work)} unique canonical SMILES")

    receptor_pdbqt = Path(args.receptor_pdbqt)
    ref_lig_sdf = Path(args.ref_lig_sdf)
    gnina_bin = Path(args.gnina)

    rows = []
    if args.jobs <= 1:
        for _, r in tqdm(df_work.iterrows(), total=len(df_work), desc="Docking"):
            out = dock_one(r, receptor_pdbqt, ref_lig_sdf, gnina_bin, lig3d_dir, gnina_out,
                           args.exhaustiveness, args.num_modes, args.seed)
            if out:
                rows.append(out)
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = []
            for _, r in df_work.iterrows():
                futs.append(ex.submit(
                    dock_one, r.to_dict(),
                    receptor_pdbqt, ref_lig_sdf, gnina_bin, lig3d_dir, gnina_out,
                    args.exhaustiveness, args.num_modes, args.seed
                ))
            for fut in tqdm(as_completed(futs), total=len(futs), desc="Docking (parallel)"):
                out = fut.result()
                if out:
                    rows.append(out)

    df_final = pd.DataFrame(rows)
    csv_path = dataset_dir/"gba1_gnina_physchem_dataset.csv"
    df_final.to_csv(csv_path, index=False)
    print("[OK] Saved:", csv_path)
    print("[INFO] Shape:", df_final.shape)

if __name__ == "__main__":
    main()
