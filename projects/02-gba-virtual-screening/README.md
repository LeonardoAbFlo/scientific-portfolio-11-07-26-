# GBA virtual-screening and affinity-modeling pipeline

This workflow curates human GBA1 activity measurements from ChEMBL, prepares a
protein structure, docks ligands with GNINA, builds a physicochemical feature
table, and trains a regression model for affinity prediction.

## Structure

- `config/config.yaml` — target, structure, docking, and split settings
- `bin/03_fetch_chembl_activities.py` — activity retrieval and unit-aware aggregation
- `bin/04_prepare_gba_receptor.py` — PDB download and receptor/reference-ligand setup
- `bin/05_build_docking_dataset.py` — ligand preparation, docking, and features
- `bin/06_train_affinity_model.py` — grouped data split, model fitting, and metrics
- `bin/07_predict_affinity.py` — prediction from a trained artifact
- `run_virtual_screening_pipeline.sh` — complete workflow

## Environment

The supplied Conda specification is in `env/ML.yml`. GNINA is downloaded at run
time by `bin/02_download_gnina.sh` and is not stored in this portfolio.

```bash
conda env create -f env/ML.yml
conda activate ML
bash run_virtual_screening_pipeline.sh
```

Generated datasets, docking poses, and model artifacts are written to `outputs/`
and intentionally excluded from the portfolio.

