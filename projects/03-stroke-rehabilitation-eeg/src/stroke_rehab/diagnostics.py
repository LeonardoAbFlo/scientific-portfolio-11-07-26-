"""Trigger, electrode-activity and time-course diagnostics."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .config import ACTIVITY_CHANNELS, FIXED_BAND, FIXED_WINDOW, RUN_TYPES, STAGES, SUBJECTS
from .epochs import find_trigger_onsets, load_preprocess_epoch
from .io import ensure_output_dirs, file_path, load_mat, save_dataframe


def resolved_channel_names(n_channels: int, preferred: Sequence[str]) -> list[str]:
    if n_channels == len(preferred):
        return list(preferred)
    return [f"Ch{index + 1}" for index in range(n_channels)]


def build_trigger_diagnostics(
    data_dir: str | Path,
    subjects: Iterable[str] = SUBJECTS,
    stages: Iterable[str] = STAGES,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for subject in subjects:
        for stage in stages:
            for run_type in RUN_TYPES:
                path = file_path(data_dir, subject, stage, run_type)
                if not path.exists():
                    continue
                fs, _, trigger = load_mat(path)
                onsets, labels = find_trigger_onsets(trigger)
                for trial_number, (onset, label) in enumerate(zip(onsets, labels), start=1):
                    rows.append(
                        {
                            "subject": subject,
                            "stage": stage.upper(),
                            "run_type": run_type.upper(),
                            "file": path.name,
                            "fs": fs,
                            "trial_number": trial_number,
                            "onset_sample": int(onset),
                            "onset_time_s": float(onset / fs),
                            "label_value": int(label),
                            "imagery_class": (
                                "Left imagery (+1)" if int(label) == 1 else "Right imagery (-1)"
                            ),
                        }
                    )
    return pd.DataFrame(rows)


def build_electrode_activity_summary(
    data_dir: str | Path,
    subjects: Iterable[str] = SUBJECTS,
    stages: Iterable[str] = STAGES,
    window: tuple[float, float] = FIXED_WINDOW,
    band: tuple[float, float] = FIXED_BAND,
    channel_names: Sequence[str] = ACTIVITY_CHANNELS,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for subject in subjects:
        for stage in stages:
            for run_type in RUN_TYPES:
                path = file_path(data_dir, subject, stage, run_type)
                if not path.exists():
                    continue
                _, epochs, labels = load_preprocess_epoch(path, window=window, band=band)
                names = resolved_channel_names(epochs.shape[2], channel_names)
                for label_value, class_name in (
                    (1, "Left imagery (+1)"),
                    (-1, "Right imagery (-1)"),
                ):
                    class_epochs = epochs[labels == label_value]
                    if class_epochs.size == 0:
                        continue
                    for channel_index, channel_name in enumerate(names):
                        values = class_epochs[:, :, channel_index]
                        rows.append(
                            {
                                "subject": subject,
                                "stage": stage.upper(),
                                "run_type": run_type.upper(),
                                "file": path.name,
                                "window": f"{window[0]}-{window[1]} s",
                                "band": f"{band[0]}-{band[1]} Hz",
                                "imagery_class": class_name,
                                "channel": channel_name,
                                "channel_index": channel_index + 1,
                                "mean_absolute_amplitude": float(np.mean(np.abs(values))),
                                "rms_amplitude": float(np.sqrt(np.mean(values**2))),
                                "mean_log_variance": float(
                                    np.mean(np.log(np.var(values, axis=1) + 1e-10))
                                ),
                                "mean_variance": float(np.mean(np.var(values, axis=1))),
                            }
                        )
    return pd.DataFrame(rows)


def build_electrode_timecourse_summary(
    data_dir: str | Path,
    subjects: Iterable[str] = SUBJECTS,
    stages: Iterable[str] = STAGES,
    window: tuple[float, float] = FIXED_WINDOW,
    band: tuple[float, float] = FIXED_BAND,
    n_time_points: int = 200,
    channel_names: Sequence[str] = ACTIVITY_CHANNELS,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for subject in subjects:
        for stage in stages:
            path = file_path(data_dir, subject, stage, "test")
            if not path.exists():
                continue
            _, epochs, labels = load_preprocess_epoch(path, window=window, band=band)
            names = resolved_channel_names(epochs.shape[2], channel_names)
            original_time = np.linspace(window[0], window[1], epochs.shape[1])
            target_time = np.linspace(window[0], window[1], n_time_points)
            for label_value, class_name in (
                (1, "Left imagery (+1)"),
                (-1, "Right imagery (-1)"),
            ):
                class_epochs = epochs[labels == label_value]
                if class_epochs.size == 0:
                    continue
                mean_timecourse = class_epochs.mean(axis=0)
                for channel_index, channel_name in enumerate(names):
                    interpolated = np.interp(
                        target_time, original_time, mean_timecourse[:, channel_index]
                    )
                    for time_value, amplitude in zip(target_time, interpolated):
                        rows.append(
                            {
                                "subject": subject,
                                "stage": stage.upper(),
                                "run_type": "TEST",
                                "imagery_class": class_name,
                                "channel": channel_name,
                                "channel_index": channel_index + 1,
                                "time_s": float(time_value),
                                "mean_amplitude": float(amplitude),
                            }
                        )
    return pd.DataFrame(rows)


def run_diagnostics(
    data_dir: str | Path,
    output_dir: str | Path,
    subjects: Iterable[str] = SUBJECTS,
    stages: Iterable[str] = STAGES,
) -> dict[str, pd.DataFrame]:
    directories = ensure_output_dirs(output_dir)
    diagnostics_dir = directories["diagnostics"]
    trigger_df = build_trigger_diagnostics(data_dir, subjects=subjects, stages=stages)
    activity_df = build_electrode_activity_summary(data_dir, subjects=subjects, stages=stages)
    timecourse_df = build_electrode_timecourse_summary(data_dir, subjects=subjects, stages=stages)

    save_dataframe(trigger_df, diagnostics_dir / "trigger_onset_diagnostics.csv")
    save_dataframe(activity_df, diagnostics_dir / "electrode_activity_summary_fixed_window.csv")
    save_dataframe(
        timecourse_df,
        diagnostics_dir / "electrode_timecourse_summary_test_fixed_window.csv",
    )
    return {
        "triggers": trigger_df,
        "electrode_activity": activity_df,
        "electrode_timecourse": timecourse_df,
    }
