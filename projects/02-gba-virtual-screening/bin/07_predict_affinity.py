\
#!/usr/bin/env python3
import argparse, json
import numpy as np
import pandas as pd
import joblib
import xgboost as xgb
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors

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

def compute_physchem(smiles, feature_cols):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    d = {}
    for k in feature_cols:
        fn = PHYS_DESC_FUNCS.get(k)
        if fn is None:
            raise ValueError(f"Unknown feature: {k}")
        d[k] = float(fn(m))
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--scaler", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--label", default="minimizedAffinity")
    ap.add_argument("--smiles", nargs="+", required=True)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    bst = xgb.Booster()
    bst.load_model(args.model)
    scaler = joblib.load(args.scaler)
    feature_cols = json.loads(open(args.features,"r").read())

    feats = []
    good = []
    for smi in args.smiles:
        smi2 = canonicalize_smiles(smi)
        if smi2 is None:
            continue
        d = compute_physchem(smi2, feature_cols)
        feats.append([d[c] for c in feature_cols])
        good.append(smi2)

    X = np.array(feats, dtype=np.float32)
    X = scaler.transform(X)
    pred = bst.predict(xgb.DMatrix(X))
    out = pd.DataFrame({"smiles": good, "pred_"+args.label: pred}).sort_values("pred_"+args.label)

    if args.out:
        out.to_csv(args.out, index=False)
        print("[OK] Saved:", args.out)
    else:
        print(out.to_string(index=False))

if __name__ == "__main__":
    main()
