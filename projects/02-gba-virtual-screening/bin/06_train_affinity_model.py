#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

import joblib
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import xgboost as xgb

def bemis_murcko_scaffold(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    scaf = MurckoScaffold.GetScaffoldForMol(m)
    if scaf is None:
        return None
    return Chem.MolToSmiles(scaf, isomericSmiles=False)

def scaffold_split(df, frac_train=0.8, frac_val=0.1, seed=0):
    df = df.copy()
    df["scaffold"] = df["smiles"].apply(bemis_murcko_scaffold)
    df = df.dropna(subset=["scaffold"]).reset_index(drop=True)

    scaffolds = df["scaffold"].value_counts().index.tolist()
    rng = np.random.default_rng(seed)
    rng.shuffle(scaffolds)

    n = len(df)
    train_cut = int(frac_train * n)
    val_cut   = int((frac_train + frac_val) * n)

    train_idx, val_idx, test_idx = [], [], []
    count = 0
    for sc in scaffolds:
        idx = df.index[df["scaffold"] == sc].tolist()
        if count < train_cut:
            train_idx += idx
        elif count < val_cut:
            val_idx += idx
        else:
            test_idx += idx
        count += len(idx)

    return df.loc[train_idx], df.loc[val_idx], df.loc[test_idx]

def train_xgb_gpu_if_possible(X_train, y_train, X_val, y_val, seed=0):
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dval   = xgb.DMatrix(X_val, label=y_val)
    evals = [(dtrain, "train"), (dval, "val")]

    base_params = {
        "objective": "reg:squarederror",
        "eval_metric": "rmse",
        "max_depth": 8,
        "eta": 0.05,
        "subsample": 0.85,
        "colsample_bytree": 0.85,
        "min_child_weight": 1.0,
        "lambda": 1.0,
        "alpha": 0.0,
        "seed": seed,
    }

    params_modern = dict(base_params)
    params_modern.update({"tree_method": "hist", "device": "cuda"})
    try:
        bst = xgb.train(params_modern, dtrain, num_boost_round=5000, evals=evals,
                        early_stopping_rounds=100, verbose_eval=50)
        return bst, "gpu(device=cuda)"
    except Exception as e1:
        print("⚠️ device=cuda failed:", str(e1).splitlines()[0])

    params_legacy = dict(base_params)
    params_legacy.update({"tree_method": "gpu_hist"})
    try:
        bst = xgb.train(params_legacy, dtrain, num_boost_round=5000, evals=evals,
                        early_stopping_rounds=100, verbose_eval=50)
        return bst, "gpu(gpu_hist)"
    except Exception as e2:
        print("⚠️ gpu_hist failed:", str(e2).splitlines()[0])

    params_cpu = dict(base_params)
    params_cpu.update({"tree_method": "hist"})
    bst = xgb.train(params_cpu, dtrain, num_boost_round=5000, evals=evals,
                    early_stopping_rounds=100, verbose_eval=50)
    return bst, "cpu"

def eval_regression(bst, X, y, name="set"):
    d = xgb.DMatrix(X)
    pred = bst.predict(d)
    mae = mean_absolute_error(y, pred)
    rmse = math.sqrt(mean_squared_error(y, pred))
    r2 = r2_score(y, pred)
    spr = spearmanr(y, pred).correlation
    print(f"[{name}] MAE={mae:.4f}  RMSE={rmse:.4f}  R2={r2:.4f}  Spearman={spr:.4f}")
    return pred

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--label", default="minimizedAffinity")
    ap.add_argument("--split-seed", type=int, default=0)
    ap.add_argument("--frac-train", type=float, default=0.8)
    ap.add_argument("--frac-val", type=float, default=0.1)
    ap.add_argument("--features", nargs="*", default=[
        "MolWt","LogP","TPSA","HBD","HBA","RotB","RingCount","AromaticRings","AliphaticRings",
        "FractionCSP3","HeavyAtomCount","NHOHCount","NOCount","MolMR"
    ])
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.dataset).dropna(subset=[args.label]).reset_index(drop=True)
    train_df, val_df, test_df = scaffold_split(df, frac_train=args.frac_train, frac_val=args.frac_val, seed=args.split_seed)
    print("[INFO] Train/Val/Test:", train_df.shape, val_df.shape, test_df.shape)

    feature_cols = args.features
    X_train = train_df[feature_cols].values.astype(np.float32)
    X_val   = val_df[feature_cols].values.astype(np.float32)
    X_test  = test_df[feature_cols].values.astype(np.float32)

    y_train = train_df[args.label].values.astype(np.float32)
    y_val   = val_df[args.label].values.astype(np.float32)
    y_test  = test_df[args.label].values.astype(np.float32)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val   = scaler.transform(X_val)
    X_test  = scaler.transform(X_test)

    bst, mode = train_xgb_gpu_if_possible(X_train, y_train, X_val, y_val, seed=args.split_seed)
    print("[INFO] Training mode:", mode)

    _ = eval_regression(bst, X_train, y_train, "train")
    _ = eval_regression(bst, X_val, y_val, "val")
    _ = eval_regression(bst, X_test, y_test, "test")

    model_path = outdir/"xgb_gba1_gnina.json"
    scaler_path = outdir/"scaler_physchem.joblib"
    feats_path = outdir/"feature_cols.json"

    bst.save_model(model_path)
    joblib.dump(scaler, scaler_path)
    feats_path.write_text(json.dumps(feature_cols, indent=2))
    print("[OK] Saved:", model_path, scaler_path, feats_path)

if __name__ == "__main__":
    main()
