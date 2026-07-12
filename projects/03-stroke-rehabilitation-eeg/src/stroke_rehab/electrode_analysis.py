"""Permutation-based channel contribution for the fixed CSP+LDA model."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score

from .config import (
    CONTRIBUTION_CHANNELS,
    FIXED_BAND,
    FIXED_N_COMPONENTS,
    FIXED_WINDOW,
    RANDOM_STATE,
    STAGES,
    SUBJECTS,
)
from .epochs import load_preprocess_epoch
from .features import fit_csp, transform_csp
from .io import ensure_output_dirs, file_path, save_dataframe
from .modeling import get_classifier


def compute_fixed_csp_lda_channel_contribution(
    train_file: str | Path,
    test_file: str | Path,
    subject: str,
    stage: str,
    channel_names: Sequence[str] = CONTRIBUTION_CHANNELS,
    n_repeats: int = 30,
    random_state: int = RANDOM_STATE,
) -> pd.DataFrame:
    """Estimate contribution as baseline accuracy minus shuffled-channel accuracy."""
    rng = np.random.default_rng(random_state)
    _, train_epochs, y_train = load_preprocess_epoch(
        train_file, window=FIXED_WINDOW, band=FIXED_BAND
    )
    _, test_epochs, y_test = load_preprocess_epoch(
        test_file, window=FIXED_WINDOW, band=FIXED_BAND
    )
    n_channels = train_epochs.shape[2]
    if n_channels != len(channel_names):
        raise ValueError(
            f"Expected {len(channel_names)} channels but found {n_channels}. "
            "Update CONTRIBUTION_CHANNELS to match the acquisition order."
        )

    weights = fit_csp(train_epochs, y_train, n_components=FIXED_N_COMPONENTS)
    x_train = transform_csp(train_epochs, weights)
    x_test = transform_csp(test_epochs, weights)
    classifier = get_classifier("lda")
    classifier.fit(x_train, y_train)
    baseline_accuracy = float(accuracy_score(y_test, classifier.predict(x_test)))

    rows: list[dict[str, object]] = []
    for channel_index, channel_name in enumerate(channel_names):
        shuffled_accuracies: list[float] = []
        for _ in range(n_repeats):
            shuffled_epochs = test_epochs.copy()
            trial_order = rng.permutation(shuffled_epochs.shape[0])
            shuffled_epochs[:, :, channel_index] = shuffled_epochs[
                trial_order, :, channel_index
            ]
            shuffled_features = transform_csp(shuffled_epochs, weights)
            shuffled_accuracies.append(
                float(accuracy_score(y_test, classifier.predict(shuffled_features)))
            )

        shuffled_mean = float(np.mean(shuffled_accuracies))
        shuffled_sd = float(np.std(shuffled_accuracies))
        contribution = baseline_accuracy * 100 - shuffled_mean * 100
        rows.append(
            {
                "subject": subject,
                "stage": stage.upper(),
                "model": "Fixed CSP+LDA",
                "channel": channel_name,
                "channel_index": channel_index + 1,
                "baseline_accuracy_percent": baseline_accuracy * 100,
                "shuffled_accuracy_percent_mean": shuffled_mean * 100,
                "shuffled_accuracy_percent_sd": shuffled_sd * 100,
                "contribution_pp": contribution,
                "contribution_pp_positive": max(contribution, 0.0),
            }
        )
    return pd.DataFrame(rows)


def compare_pre_post(contributions: pd.DataFrame) -> pd.DataFrame:
    pre = contributions[contributions["stage"] == "PRE"].copy()
    post = contributions[contributions["stage"] == "POST"].copy()
    comparison = pre.merge(
        post,
        on=["subject", "channel", "channel_index", "model"],
        suffixes=("_pre", "_post"),
    )
    comparison["contribution_change_pp"] = (
        comparison["contribution_pp_post"] - comparison["contribution_pp_pre"]
    )
    comparison["contribution_change_norm"] = (
        comparison["contribution_norm_post"] - comparison["contribution_norm_pre"]
    )
    comparison["change_direction"] = np.select(
        [
            comparison["contribution_change_pp"] > 0,
            comparison["contribution_change_pp"] < 0,
        ],
        ["increased", "decreased"],
        default="stable",
    )
    return comparison


def run_channel_contributions(
    data_dir: str | Path,
    output_dir: str | Path,
    subjects: Iterable[str] = SUBJECTS,
    stages: Iterable[str] = STAGES,
    n_repeats: int = 30,
) -> dict[str, pd.DataFrame]:
    directories = ensure_output_dirs(output_dir)
    frames: list[pd.DataFrame] = []
    for subject in subjects:
        for stage in stages:
            train_file = file_path(data_dir, subject, stage, "training")
            test_file = file_path(data_dir, subject, stage, "test")
            if not train_file.exists() or not test_file.exists():
                continue
            print(f"Channel contributions: {subject} {stage.upper()}")
            frames.append(
                compute_fixed_csp_lda_channel_contribution(
                    train_file,
                    test_file,
                    subject=subject,
                    stage=stage,
                    n_repeats=n_repeats,
                )
            )
    if not frames:
        raise RuntimeError("No complete training/test pairs available for channel analysis.")

    contributions = pd.concat(frames, ignore_index=True)
    contributions["contribution_norm"] = contributions.groupby(
        ["subject", "stage"]
    )["contribution_pp_positive"].transform(
        lambda values: values / values.max() if values.max() > 0 else 0.0
    )
    comparison = compare_pre_post(contributions)

    output = directories["electrodes"]
    save_dataframe(
        contributions,
        output / "fixed_csp_lda_channel_contributions_long.csv",
        output / "fixed_csp_lda_channel_contributions_long.json",
    )
    save_dataframe(
        comparison,
        output / "fixed_csp_lda_channel_contributions_pre_post_compare.csv",
        output / "fixed_csp_lda_channel_contributions_pre_post_compare.json",
    )
    return {"contributions": contributions, "comparison": comparison}
