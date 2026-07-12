#!/usr/bin/env bash
set -euo pipefail

CONFIG="config/config.yaml"

get_cfg () {
  local key="$1"
  python - "$key" <<'PY'
import sys, yaml
from pathlib import Path

key = sys.argv[1]
cfg_path = Path("config/config.yaml")
cfg = yaml.safe_load(cfg_path.read_text())

val = cfg.get(key)
if isinstance(val, str):
    print(val)
elif val is None:
    print("")
else:
    print(val)
PY
}

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ML
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

python - <<'PY' >/dev/null 2>&1 || pip -q install pyyaml
import yaml
PY

if [[ ! -f "$CONFIG" ]]; then
  echo "[ERR] Missing config file: $CONFIG"
  exit 1
fi

bash bin/01_check_environment.sh
bash bin/02_download_gnina.sh

TARGET="$(get_cfg target_chembl_id)"
PDB_ID="$(get_cfg pdb_id)"
N_LIGANDS="$(get_cfg n_ligands)"
EXH="$(get_cfg exhaustiveness)"
NM="$(get_cfg num_modes)"
SEED="$(get_cfg seed)"
JOBS="$(get_cfg jobs)"
LABEL="$(get_cfg label)"
SPLIT_SEED="$(get_cfg split_seed)"
FRAC_TRAIN="$(get_cfg frac_train)"
FRAC_VAL="$(get_cfg frac_val)"

mkdir -p outputs/chembl outputs/receptor outputs/datasets outputs/model outputs/docking

echo "[STEP] Fetch ChEMBL activities for ${TARGET}"
python bin/03_fetch_chembl_activities.py --target "${TARGET}" --outdir outputs/chembl

echo "[STEP] Prepare receptor from PDB ${PDB_ID}"
python bin/04_prepare_gba_receptor.py --pdb-id "${PDB_ID}" --outdir outputs/receptor

GNINA_BIN="tools/gnina/gnina"
RECEPTOR_PDBQT="outputs/receptor/${PDB_ID}_protein.pdbqt"
REF_LIG_SDF="outputs/receptor/${PDB_ID}_ref_lig.sdf"

if [[ ! -x "$GNINA_BIN" ]]; then
  echo "[ERR] GNINA binary not executable or missing: $GNINA_BIN"
  echo "      Try: bash bin/02_download_gnina.sh"
  exit 1
fi
if [[ ! -f "$RECEPTOR_PDBQT" ]]; then
  echo "[ERR] Missing receptor pdbqt: $RECEPTOR_PDBQT"
  exit 1
fi
if [[ ! -f "$REF_LIG_SDF" ]]; then
  echo "[ERR] Missing ref ligand sdf: $REF_LIG_SDF"
  exit 1
fi

echo "[STEP] Dock + build dataset"
python bin/05_build_docking_dataset.py \
  --agg-csv outputs/chembl/gba1_agg_medians.csv \
  --receptor-pdbqt "$RECEPTOR_PDBQT" \
  --ref-lig-sdf "$REF_LIG_SDF" \
  --gnina "$GNINA_BIN" \
  --outdir outputs \
  --n-ligands "$N_LIGANDS" \
  --exhaustiveness "$EXH" \
  --num-modes "$NM" \
  --seed "$SEED" \
  --jobs "$JOBS"

echo "[STEP] Train model"
python bin/06_train_affinity_model.py \
  --dataset outputs/datasets/gba1_gnina_physchem_dataset.csv \
  --outdir outputs/model \
  --label "$LABEL" \
  --split-seed "$SPLIT_SEED" \
  --frac-train "$FRAC_TRAIN" \
  --frac-val "$FRAC_VAL"

echo "[OK] Pipeline completed."
