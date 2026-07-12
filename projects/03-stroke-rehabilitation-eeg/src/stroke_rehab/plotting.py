"""Matplotlib figures for the main decoding outputs."""

from __future__ import annotations

from ast import literal_eval
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import ConfusionMatrixDisplay

from .config import LABEL_NAMES, STAGES, SUBJECTS
from .io import ensure_output_dirs


def _read_results(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    if "confusion_matrix" in frame.columns:
        frame["confusion_matrix"] = frame["confusion_matrix"].apply(literal_eval)
    return frame


def plot_accuracy_lines(
    results: pd.DataFrame,
    output_path: Path,
    title: str,
    subjects: Iterable[str] = SUBJECTS,
) -> None:
    fig, axis = plt.subplots(figsize=(8, 5))
    for subject in subjects:
        subset = results[results["subject"] == subject].copy()
        if subset.empty:
            continue
        subset["stage"] = pd.Categorical(subset["stage"], categories=STAGES, ordered=True)
        subset = subset.sort_values("stage")
        axis.plot(subset["stage"].astype(str), subset["test_accuracy_percent"], marker="o", label=subject)
    axis.axhline(50, linestyle="--", linewidth=1, label="Chance level (50%)")
    axis.set_ylim(0, 105)
    axis.set_ylabel("Held-out test accuracy (%)")
    axis.set_xlabel("Session")
    axis.set_title(title)
    axis.legend()
    axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_change_bars(change: pd.DataFrame, output_path: Path, title: str) -> None:
    fig, axis = plt.subplots(figsize=(7, 5))
    axis.bar(change["subject"], change["post_minus_pre_percent"])
    axis.axhline(0, linewidth=1)
    axis.set_ylabel("POST - PRE accuracy (percentage points)")
    axis.set_xlabel("Subject")
    axis.set_title(title)
    axis.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_confusion_matrices(
    results: pd.DataFrame,
    output_dir: Path,
    model_prefix: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for _, row in results.iterrows():
        matrix = np.asarray(row["confusion_matrix"])
        figure, axis = plt.subplots(figsize=(7, 5))
        ConfusionMatrixDisplay(matrix, display_labels=LABEL_NAMES).plot(
            ax=axis, cmap="Blues", values_format="d", colorbar=True
        )
        axis.set_title(
            f"{row['subject']} {str(row['stage']).upper()} — "
            f"{model_prefix} — {row['test_accuracy_percent']:.1f}%"
        )
        figure.tight_layout()
        filename = f"{model_prefix.lower().replace(' ', '_').replace('+', '_')}_{row['subject']}_{row['stage']}.png"
        figure.savefig(output_dir / filename, dpi=300)
        plt.close(figure)


def make_all_python_figures(output_dir: str | Path) -> None:
    directories = ensure_output_dirs(output_dir)
    tables = directories["tables"]
    figures = directories["figures"]

    fixed = _read_results(tables / "stroke_bci_fixed_csp_lda_results.csv")
    selected = _read_results(tables / "stroke_bci_training_cv_selected_results.csv")
    fixed_change = pd.read_csv(tables / "stroke_bci_fixed_change_pre_post.csv")
    selected_change = pd.read_csv(tables / "stroke_bci_selected_change_pre_post.csv")

    plot_accuracy_lines(
        fixed,
        figures / "fixed_csp_lda_pre_post_accuracy.png",
        "Fixed CSP+LDA motor imagery decoding: PRE vs POST",
    )
    plot_accuracy_lines(
        selected,
        figures / "training_cv_selected_pre_post_accuracy.png",
        "Training-CV selected model: PRE vs POST",
    )
    plot_change_bars(
        fixed_change,
        figures / "fixed_csp_lda_post_minus_pre.png",
        "Fixed CSP+LDA accuracy change after BCI rehabilitation",
    )
    plot_change_bars(
        selected_change,
        figures / "training_cv_selected_post_minus_pre.png",
        "Training-CV selected model change after BCI rehabilitation",
    )
    plot_confusion_matrices(fixed, figures / "confusion_matrices", "Fixed CSP+LDA")
    plot_confusion_matrices(selected, figures / "confusion_matrices", "Selected model")
