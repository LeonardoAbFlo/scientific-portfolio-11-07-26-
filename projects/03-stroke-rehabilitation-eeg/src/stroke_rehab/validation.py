"""Input dataset validation."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

from .config import RUN_TYPES, STAGES, SUBJECTS
from .epochs import find_trigger_onsets
from .io import ensure_output_dirs, file_path, load_mat, save_dataframe


def validate_dataset(
    data_dir: str | Path,
    output_dir: str | Path,
    subjects: Iterable[str] = SUBJECTS,
    stages: Iterable[str] = STAGES,
    strict: bool = True,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    missing: list[str] = []
    for subject in subjects:
        for stage in stages:
            for run_type in RUN_TYPES:
                path = file_path(data_dir, subject, stage, run_type)
                if not path.exists():
                    missing.append(path.name)
                    rows.append(
                        {
                            "subject": subject,
                            "stage": stage,
                            "run_type": run_type,
                            "file": path.name,
                            "status": "missing",
                        }
                    )
                    continue
                try:
                    fs, eeg, trigger = load_mat(path)
                    onsets, labels = find_trigger_onsets(trigger)
                    rows.append(
                        {
                            "subject": subject,
                            "stage": stage,
                            "run_type": run_type,
                            "file": path.name,
                            "status": "ok",
                            "fs_hz": fs,
                            "n_samples": eeg.shape[0],
                            "n_channels": eeg.shape[1],
                            "n_trials": len(onsets),
                            "n_left": int((labels == 1).sum()),
                            "n_right": int((labels == -1).sum()),
                        }
                    )
                except Exception as exc:
                    rows.append(
                        {
                            "subject": subject,
                            "stage": stage,
                            "run_type": run_type,
                            "file": path.name,
                            "status": "invalid",
                            "error": str(exc),
                        }
                    )

    report = pd.DataFrame(rows)
    diagnostics = ensure_output_dirs(output_dir)["diagnostics"]
    save_dataframe(report, diagnostics / "data_validation_report.csv")
    invalid = report[report["status"] != "ok"]
    if strict and not invalid.empty:
        raise RuntimeError(
            "Dataset validation failed. See results/diagnostics/data_validation_report.csv."
        )
    return report
