"""End-to-end decoding evaluation and result export."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .config import FIXED_N_COMPONENTS, FIXED_WINDOW, LABEL_NAMES, STAGES, SUBJECTS
from .io import ensure_output_dirs, file_path, save_dataframe
from .modeling import evaluate_configuration, evaluate_fixed_csp_lda, select_training_only_configuration


def build_change_table(
    results_df: pd.DataFrame,
    subjects: Iterable[str] = SUBJECTS,
    accuracy_col: str = "test_accuracy_percent",
) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    for subject in subjects:
        subset = results_df[results_df["subject"] == subject]
        if {"pre", "post"}.issubset(set(subset["stage"])):
            pre = float(subset.loc[subset["stage"] == "pre", accuracy_col].iloc[0])
            post = float(subset.loc[subset["stage"] == "post", accuracy_col].iloc[0])
            rows.append(
                {
                    "subject": subject,
                    "pre_accuracy_percent": pre,
                    "post_accuracy_percent": post,
                    "post_minus_pre_percent": post - pre,
                }
            )
    return pd.DataFrame(rows)


def export_confusion_matrix_tables(
    results_df: pd.DataFrame,
    model_family: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in results_df.iterrows():
        matrix = np.asarray(row["confusion_matrix"])
        for true_index, true_label in enumerate(LABEL_NAMES):
            for predicted_index, predicted_label in enumerate(LABEL_NAMES):
                rows.append(
                    {
                        "subject": row["subject"],
                        "stage": str(row["stage"]).upper(),
                        "model_family": model_family,
                        "method": row.get("method", row.get("selected_method", "")),
                        "window": str(row.get("window", row.get("selected_window", ""))),
                        "test_accuracy_percent": float(row["test_accuracy_percent"]),
                        "true_label": true_label,
                        "predicted_label": predicted_label,
                        "n": int(matrix[true_index, predicted_index]),
                    }
                )
    return pd.DataFrame(rows)


def run_decoding(
    data_dir: str | Path,
    output_dir: str | Path,
    subjects: Iterable[str] = SUBJECTS,
    stages: Iterable[str] = STAGES,
) -> dict[str, pd.DataFrame]:
    """Run training-only configuration selection and held-out test evaluation."""
    directories = ensure_output_dirs(output_dir)
    selected_rows: list[dict[str, object]] = []
    fixed_rows: list[dict[str, object]] = []
    configuration_tables: list[pd.DataFrame] = []

    for subject in subjects:
        for stage in stages:
            train_file = file_path(data_dir, subject, stage, "training")
            test_file = file_path(data_dir, subject, stage, "test")
            if not train_file.exists() or not test_file.exists():
                print(f"Skipping missing pair: {train_file.name}, {test_file.name}")
                continue

            print(f"Evaluating {subject} {stage.upper()}")
            best_config, grid = select_training_only_configuration(train_file)
            grid["subject"] = subject
            grid["stage"] = stage
            configuration_tables.append(grid)

            selected = evaluate_configuration(train_file, test_file, best_config)
            fixed = evaluate_fixed_csp_lda(train_file, test_file)

            selected_rows.append(
                {
                    "subject": subject,
                    "stage": stage,
                    "selected_window": str(best_config["window"]),
                    "selected_method": best_config["method"],
                    "selected_classifier": best_config["classifier"],
                    "selected_n_components": best_config["n_components"],
                    "training_cv_accuracy_mean": best_config["cv_accuracy_mean"],
                    "training_cv_accuracy_sd": best_config["cv_accuracy_sd"],
                    "test_accuracy": selected["accuracy"],
                    "test_accuracy_percent": selected["accuracy_percent"],
                    "n_test_trials": selected["n_test_trials"],
                    "n_correct": selected["n_correct"],
                    "p_vs_chance_0.5": selected["p_vs_chance_0.5"],
                    "confusion_matrix": selected["confusion_matrix"],
                }
            )
            fixed_rows.append(
                {
                    "subject": subject,
                    "stage": stage,
                    "method": "Fixed CSP+LDA",
                    "window": str(FIXED_WINDOW),
                    "n_components": FIXED_N_COMPONENTS,
                    "test_accuracy": fixed["accuracy"],
                    "test_accuracy_percent": fixed["accuracy_percent"],
                    "n_test_trials": fixed["n_test_trials"],
                    "n_correct": fixed["n_correct"],
                    "p_vs_chance_0.5": fixed["p_vs_chance_0.5"],
                    "confusion_matrix": fixed["confusion_matrix"],
                }
            )
            print(
                f"  selected={best_config['method']} {best_config['classifier']} "
                f"window={best_config['window']} | test={selected['accuracy_percent']:.1f}%"
            )
            print(f"  fixed CSP+LDA | test={fixed['accuracy_percent']:.1f}%")

    selected_df = pd.DataFrame(selected_rows)
    fixed_df = pd.DataFrame(fixed_rows)
    grid_df = (
        pd.concat(configuration_tables, ignore_index=True)
        if configuration_tables
        else pd.DataFrame()
    )
    if selected_df.empty or fixed_df.empty:
        raise RuntimeError("No complete training/test pairs were evaluated.")

    selected_change = build_change_table(selected_df, subjects=subjects)
    fixed_change = build_change_table(fixed_df, subjects=subjects)
    selected_cm = export_confusion_matrix_tables(
        selected_df, "Training-CV selected CSP/FBCSP model"
    )
    fixed_cm = export_confusion_matrix_tables(fixed_df, "Fixed CSP+LDA")

    tables = directories["tables"]
    save_dataframe(
        selected_df,
        tables / "stroke_bci_training_cv_selected_results.csv",
        tables / "stroke_bci_training_cv_selected_results.json",
    )
    save_dataframe(
        fixed_df,
        tables / "stroke_bci_fixed_csp_lda_results.csv",
        tables / "stroke_bci_fixed_csp_lda_results.json",
    )
    save_dataframe(grid_df, tables / "stroke_bci_training_cv_grid.csv")
    save_dataframe(selected_change, tables / "stroke_bci_selected_change_pre_post.csv")
    save_dataframe(fixed_change, tables / "stroke_bci_fixed_change_pre_post.csv")
    save_dataframe(
        selected_cm, tables / "training_cv_selected_confusion_matrices_long.csv"
    )
    save_dataframe(fixed_cm, tables / "fixed_csp_lda_confusion_matrices_long.csv")

    return {
        "selected": selected_df,
        "fixed": fixed_df,
        "grid": grid_df,
        "selected_change": selected_change,
        "fixed_change": fixed_change,
    }
