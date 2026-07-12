#!/usr/bin/env bash
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ML

echo "[OK] Conda env: $(python -c 'import sys; print(sys.executable)')"
python -c "import numpy, pandas, scipy, sklearn, tqdm, requests, Bio, joblib; print('[OK] python deps import')"
python -c "from rdkit import Chem; print('[OK] rdkit', Chem.__version__ if hasattr(Chem,'__version__') else '')"
python -c "import xgboost as xgb; print('[OK] xgboost', xgb.__version__)"

command -v obabel >/dev/null 2>&1 && echo "[OK] obabel found: $(command -v obabel)" || { echo "[ERR] obabel not found. Install openbabel from conda-forge."; exit 1; }

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi || true
else
  echo "[INFO] nvidia-smi not found (CPU-only is fine)."
fi
