# Stroke Rehabilitation EEG Decoding Pipeline

A script-based, reproducible version of the original Jupyter notebook for binary left-versus-right motor-imagery decoding before and after stroke rehabilitation.

## What the pipeline does

1. Validates the expected MATLAB EEG files and trigger labels.
2. Removes the DC component, applies a 50 Hz notch filter, and band-pass filters EEG data.
3. Detects `+1` and `-1` trigger onsets and extracts motor-imagery epochs.
4. Compares log-variance, CSP, and filter-bank CSP features with LDA and RBF-SVM classifiers.
5. Selects the model using **training data only**, then evaluates it once on the held-out test file.
6. Runs a fixed CSP(4)+LDA analysis for consistent PRE/POST comparison.
7. Exports accuracy tables, binomial tests against 50% chance, confusion matrices, trigger diagnostics, electrode summaries, and optional channel-contribution estimates.
8. Produces Python figures and optional publication-style R figures.

## Repository structure

```text
stroke-rehab-eeg-pipeline/
├── data/raw/                     # Place .mat files here; ignored by Git
├── results/                      # Generated outputs; ignored by Git
├── scripts/
│   ├── stage_01_validate_eeg_data.py
│   ├── stage_02_run_decoding_analysis.py
│   ├── stage_03_export_diagnostics.py
│   ├── stage_04_analyze_channel_contributions.py
│   ├── stage_05_create_figures.py
│   └── run_eeg_pipeline.py
├── src/stroke_rehab/             # Reusable analysis package
├── r/create_publication_figures.r       # Optional ggplot2 figures
├── tests/test_eeg_pipeline.py
├── pyproject.toml
└── requirements.txt
```

## Expected data names

The default configuration expects 12 files:

```text
P1_pre_training.mat    P1_pre_test.mat
P1_post_training.mat   P1_post_test.mat
P2_pre_training.mat    P2_pre_test.mat
P2_post_training.mat   P2_post_test.mat
P3_pre_training.mat    P3_pre_test.mat
P3_post_training.mat   P3_post_test.mat
```

Each file must contain:

- `fs`: sampling frequency
- `y`: EEG array, interpreted as samples × channels
- `trig`: sample-level trigger vector containing `+1`, `-1`, and background values

## Installation

```bash
python -m venv .venv
source .venv/bin/activate       # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -e .
```

For tests:

```bash
pip install -r requirements-dev.txt
pytest
```

## Run the complete pipeline

```bash
python scripts/run_eeg_pipeline.py \
  --data-dir data/raw \
  --output-dir results
```

The permutation-based channel analysis is slower. Skip it during development with:

```bash
python scripts/run_eeg_pipeline.py \
  --data-dir data/raw \
  --output-dir results \
  --skip-channel-contributions
```

## Run individual steps

```bash
python scripts/stage_01_validate_eeg_data.py --data-dir data/raw --output-dir results
python scripts/stage_02_run_decoding_analysis.py --data-dir data/raw --output-dir results
python scripts/stage_03_export_diagnostics.py --data-dir data/raw --output-dir results
python scripts/stage_04_analyze_channel_contributions.py --data-dir data/raw --output-dir results --repeats 30
python scripts/stage_05_create_figures.py --output-dir results
```

## Optional R figures

After the Python tables and diagnostics have been generated:

```bash
Rscript r/create_publication_figures.r results
```

The R script installs missing plotting packages and writes PDF and 600-dpi PNG figures to `results/publication_figures/`.

## Important methodological improvement

The notebook calculated CSP/FBCSP features before cross-validation. Because CSP is fitted using class labels, that can leak information from validation folds into the feature transform. This repository fits CSP separately inside every training fold, then transforms the corresponding validation fold. Final model selection still uses only the training file, and the selected model is evaluated once on the held-out test file.

## Channel-order note

The original notebook contained two different 16-channel name lists: one for activity summaries and another for the permutation-contribution analysis. Both are preserved in `src/stroke_rehab/config.py`. Confirm the true acquisition order before interpreting electrode-specific results.

## GitHub portfolio checklist

- Do not commit participant `.mat` files unless sharing is ethically and legally permitted.
- Add a short project screenshot or example result generated from de-identified/public data.
- Describe the cohort, acquisition system, ethics approval, and rehabilitation protocol in the repository only when disclosure is permitted.
- Add a license only after choosing terms compatible with the data and institutional policy.
