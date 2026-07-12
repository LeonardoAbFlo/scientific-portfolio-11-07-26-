\
#!/usr/bin/env python3
import argparse, json, time
from pathlib import Path
from functools import lru_cache

import pandas as pd
import requests
from tqdm import tqdm

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

def chembl_get(url, params=None, timeout=60, retries=5, backoff=2.0):
    headers = {"Accept":"application/json"}
    last = None
    for i in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last = e
            sleep = backoff ** i
            time.sleep(min(sleep, 30))
    raise RuntimeError(f"ChEMBL request failed after {retries} retries: {last}")

def fetch_activities_for_target(target_chembl_id, limit=200, max_pages=50):
    out = []
    offset = 0
    for _ in range(max_pages):
        data = chembl_get(
            f"{CHEMBL_BASE}/activity.json",
            params={"target_chembl_id": target_chembl_id, "limit": limit, "offset": offset},
        )
        acts = data.get("activities", [])
        if not acts:
            break
        out.extend(acts)
        page_meta = data.get("page_meta", {})
        if page_meta.get("next") is None:
            break
        offset += limit
    return out

def clean_activity_rows(acts, keep_types, keep_relations):
    rows = []
    for a in acts:
        stype = a.get("standard_type")
        units = a.get("standard_units")
        rel = a.get("standard_relation")
        val = a.get("standard_value")
        mid = a.get("molecule_chembl_id")
        assay = a.get("assay_chembl_id")

        if mid is None or stype is None or units is None or val is None:
            continue
        if stype not in keep_types:
            continue
        if rel not in keep_relations:
            continue
        try:
            val = float(val)
        except Exception:
            continue

        rows.append({
            "molecule_chembl_id": mid,
            "assay_chembl_id": assay,
            "standard_type": stype,
            "standard_relation": rel,
            "standard_units": units,
            "standard_value_nM": val,
            "pchembl_value": a.get("pchembl_value"),
        })
    return pd.DataFrame(rows)

@lru_cache(maxsize=50000)
def fetch_smiles_for_molecule(molecule_chembl_id):
    data = chembl_get(f"{CHEMBL_BASE}/molecule/{molecule_chembl_id}.json")
    ms = data.get("molecule_structures", {})
    return ms.get("canonical_smiles")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, help="ChEMBL target id (e.g., CHEMBL2179)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--limit", type=int, default=200)
    ap.add_argument("--max-pages", type=int, default=50)
    ap.add_argument("--keep-types", nargs="+", default=["IC50","Ki","Kd"])
    ap.add_argument("--keep-relations", nargs="+", default=["=","~"])
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    acts = fetch_activities_for_target(args.target, limit=args.limit, max_pages=args.max_pages)
    print(f"[INFO] Activities fetched: {len(acts)}")

    df_act = clean_activity_rows(acts, set(args.keep_types), set(args.keep_relations))
    print(f"[INFO] Filtered activity rows: {df_act.shape}")

    df_act["pchembl_value"] = pd.to_numeric(df_act["pchembl_value"], errors="coerce")

    mol_ids = sorted(df_act["molecule_chembl_id"].unique().tolist())
    print(f"[INFO] Unique molecules: {len(mol_ids)}")

    smiles_map = {}
    for mid in tqdm(mol_ids, desc="Fetch SMILES"):
        smi = fetch_smiles_for_molecule(mid)
        if smi:
            smiles_map[mid] = smi

    df_smiles = pd.DataFrame({
        "molecule_chembl_id": list(smiles_map.keys()),
        "canonical_smiles": list(smiles_map.values())
    })
    print(f"[INFO] Molecules with SMILES: {df_smiles.shape}")

    df = df_act.merge(df_smiles, on="molecule_chembl_id", how="inner")

    agg = (df.groupby(["molecule_chembl_id","canonical_smiles"], as_index=False)
             .agg(standard_value_nM_median=("standard_value_nM","median"),
                  n_activity_points=("standard_value_nM","count"),
                  pchembl_median=("pchembl_value","median")))
    agg = agg.sort_values("n_activity_points", ascending=False).reset_index(drop=True)

    df_act.to_csv(outdir/"gba1_activities_clean.csv", index=False)
    agg.to_csv(outdir/"gba1_agg_medians.csv", index=False)
    print("[OK] Saved:",
          outdir/"gba1_activities_clean.csv",
          outdir/"gba1_agg_medians.csv")

if __name__ == "__main__":
    main()
