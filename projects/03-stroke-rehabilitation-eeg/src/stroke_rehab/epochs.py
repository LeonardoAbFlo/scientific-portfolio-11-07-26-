"""Trigger detection and epoch extraction."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .io import load_mat
from .preprocessing import preprocess_continuous


def find_trigger_onsets(trig: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Find transitions into valid +1/-1 trigger blocks."""
    trig = np.asarray(trig).squeeze()
    valid = np.isin(trig, [-1, 1])
    indices = np.where(valid)[0]
    if indices.size == 0:
        raise ValueError("No valid triggers found; expected +1 and -1.")

    onsets: list[int] = []
    labels: list[int] = []
    previous = -10**12
    for index in indices:
        is_onset = index == 0 or not valid[index - 1]
        if is_onset and index - previous > 1:
            onsets.append(int(index))
            labels.append(int(trig[index]))
            previous = int(index)
    return np.asarray(onsets, dtype=int), np.asarray(labels, dtype=int)


def extract_epochs_from_continuous(
    y: np.ndarray,
    trig: np.ndarray,
    fs: int,
    window: tuple[float, float] = (2.0, 8.0),
) -> tuple[np.ndarray, np.ndarray]:
    """Extract complete epochs after each trigger onset."""
    tmin, tmax = window
    if tmax <= tmin:
        raise ValueError(f"Invalid epoch window {window}.")

    start_offset = int(round(tmin * fs))
    end_offset = int(round(tmax * fs))
    onsets, raw_labels = find_trigger_onsets(trig)

    epochs: list[np.ndarray] = []
    labels: list[int] = []
    for onset, label in zip(onsets, raw_labels):
        start = int(onset + start_offset)
        end = int(onset + end_offset)
        if start >= 0 and end <= y.shape[0]:
            epochs.append(y[start:end, :])
            labels.append(int(label))

    epoch_array = np.asarray(epochs, dtype=float)
    label_array = np.asarray(labels, dtype=int)
    if epoch_array.ndim != 3 or epoch_array.shape[0] == 0:
        raise ValueError(
            "Epoch extraction failed; expected a non-empty trials x samples x channels array."
        )
    return epoch_array, label_array


def load_preprocess_epoch(
    path: str | Path,
    window: tuple[float, float] = (2.0, 8.0),
    band: tuple[float, float] = (8.0, 30.0),
) -> tuple[int, np.ndarray, np.ndarray]:
    fs, y, trig = load_mat(path)
    filtered = preprocess_continuous(y, fs, band=band, use_notch=True)
    epochs, labels = extract_epochs_from_continuous(filtered, trig, fs, window=window)
    return fs, epochs, labels
