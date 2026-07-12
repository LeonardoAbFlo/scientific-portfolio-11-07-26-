"""Input/output helpers for MATLAB EEG files and tabular results."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import scipy.io as sio


def load_mat(path: str | Path) -> tuple[int, np.ndarray, np.ndarray]:
    """Load ``fs``, ``y`` and ``trig`` from MATLAB v7 or v7.3 files.

    Returns
    -------
    fs
        Sampling frequency in Hz.
    y
        EEG matrix with shape ``samples x channels``.
    trig
        Trigger vector with one value per sample.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    try:
        mat = sio.loadmat(path)
        fs = int(np.squeeze(mat["fs"]))
        y = np.asarray(mat["y"], dtype=float)
        trig = np.asarray(mat["trig"]).squeeze()
    except (NotImplementedError, ValueError):
        with h5py.File(path, "r") as handle:
            fs = int(np.asarray(handle["fs"]).squeeze())
            y = np.asarray(handle["y"], dtype=float)
            trig = np.asarray(handle["trig"]).squeeze()

    if y.ndim != 2:
        raise ValueError(f"Expected a 2D EEG matrix in {path}, found shape {y.shape}.")
    if y.shape[0] < y.shape[1]:
        y = y.T

    trig = np.asarray(trig).reshape(-1)
    if trig.size != y.shape[0]:
        raise ValueError(
            f"Trigger length {trig.size} does not match EEG samples {y.shape[0]} in {path}."
        )
    if fs <= 0:
        raise ValueError(f"Sampling frequency must be positive in {path}; found {fs}.")

    return fs, y, trig


def file_path(data_dir: str | Path, subject: str, stage: str, run_type: str) -> Path:
    """Construct the expected input filename."""
    return Path(data_dir) / f"{subject}_{stage}_{run_type}.mat"


def ensure_output_dirs(output_dir: str | Path) -> dict[str, Path]:
    """Create and return the standard output directories."""
    root = Path(output_dir)
    directories = {
        "root": root,
        "tables": root / "tables",
        "figures": root / "figures",
        "diagnostics": root / "diagnostics",
        "electrodes": root / "electrodes",
        "publication_figures": root / "publication_figures",
    }
    for path in directories.values():
        path.mkdir(parents=True, exist_ok=True)
    return directories


def dataframe_records(df: pd.DataFrame) -> list[dict[str, Any]]:
    """Convert a dataframe to JSON-safe records."""
    records = df.replace({np.nan: None}).to_dict(orient="records")
    for record in records:
        for key, value in list(record.items()):
            if isinstance(value, np.generic):
                record[key] = value.item()
            elif isinstance(value, np.ndarray):
                record[key] = value.tolist()
    return records


def save_dataframe(df: pd.DataFrame, csv_path: str | Path, json_path: str | Path | None = None) -> None:
    """Save a dataframe as CSV and optionally as record-oriented JSON."""
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    if json_path is not None:
        json_path = Path(json_path)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(dataframe_records(df), indent=2), encoding="utf-8")
