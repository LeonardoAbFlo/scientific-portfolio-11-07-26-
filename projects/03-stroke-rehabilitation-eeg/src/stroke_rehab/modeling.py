"""Leakage-safe model selection and held-out evaluation."""

from __future__ import annotations

from ast import literal_eval
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from .config import (
    CLASSIFIER_GRID,
    FIXED_N_COMPONENTS,
    FIXED_WINDOW,
    LABEL_ORDER,
    METHOD_GRID,
    N_COMPONENTS_GRID,
    RANDOM_STATE,
    WINDOW_GRID,
)
from .features import fit_features, load_epoch_blocks, transform_features


def get_classifier(name: str) -> Pipeline:
    if name == "lda":
        return Pipeline(
            [("scaler", StandardScaler()), ("clf", LinearDiscriminantAnalysis())]
        )
    if name == "svm":
        return Pipeline(
            [("scaler", StandardScaler()), ("clf", SVC(kernel="rbf", C=1.0, gamma="scale"))]
        )
    raise ValueError(f"Unknown classifier: {name}.")


def _minimum_class_count(labels: np.ndarray) -> int:
    _, counts = np.unique(labels, return_counts=True)
    return int(counts.min()) if counts.size else 0


def training_cv_score(
    train_file: str | Path,
    method: str = "csp",
    classifier: str = "lda",
    window: tuple[float, float] = (2.0, 8.0),
    n_components: int = 4,
    n_splits: int = 5,
) -> tuple[float, float]:
    """Score a configuration using fold-wise CSP fitting to prevent leakage."""
    epoch_blocks, labels = load_epoch_blocks(train_file, method=method, window=window)
    n_splits = min(int(n_splits), _minimum_class_count(labels))
    if n_splits < 2:
        return float("nan"), float("nan")

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_STATE)
    scores: list[float] = []
    dummy = np.zeros(labels.shape[0])
    for train_index, validation_index in cv.split(dummy, labels):
        train_blocks = [block[train_index] for block in epoch_blocks]
        validation_blocks = [block[validation_index] for block in epoch_blocks]
        y_train = labels[train_index]
        y_validation = labels[validation_index]

        x_train, weights = fit_features(
            train_blocks, y_train, method=method, n_components=n_components
        )
        x_validation = transform_features(validation_blocks, method=method, weights=weights)
        model = get_classifier(classifier)
        model.fit(x_train, y_train)
        scores.append(float(model.score(x_validation, y_validation)))

    return float(np.mean(scores)), float(np.std(scores))


def select_training_only_configuration(train_file: str | Path) -> tuple[dict[str, Any], pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    simplicity = {"logvar": 0, "csp": 1, "fbcsp": 2}

    for window in WINDOW_GRID:
        for method in METHOD_GRID:
            for classifier in CLASSIFIER_GRID:
                component_grid: tuple[int | None, ...] = (
                    (None,) if method == "logvar" else tuple(N_COMPONENTS_GRID)
                )
                for n_components in component_grid:
                    row: dict[str, Any] = {
                        "window": str(tuple(window)),
                        "method": method,
                        "classifier": classifier,
                        "n_components": np.nan if n_components is None else int(n_components),
                        "simplicity": simplicity[method],
                        "n_comp_sort": 0 if n_components is None else int(n_components),
                    }
                    try:
                        mean_accuracy, sd_accuracy = training_cv_score(
                            train_file,
                            method=method,
                            classifier=classifier,
                            window=tuple(window),
                            n_components=(
                                FIXED_N_COMPONENTS if n_components is None else int(n_components)
                            ),
                            n_splits=5,
                        )
                        row.update(
                            cv_accuracy_mean=mean_accuracy,
                            cv_accuracy_sd=sd_accuracy,
                            error=None,
                        )
                    except Exception as exc:  # grid search should record failures, not abort
                        row.update(
                            cv_accuracy_mean=np.nan,
                            cv_accuracy_sd=np.nan,
                            error=str(exc),
                        )
                    rows.append(row)

    table = pd.DataFrame(rows)
    valid = table.dropna(subset=["cv_accuracy_mean"]).copy()
    if valid.empty:
        raise ValueError("No valid training-only configuration was found.")
    valid = valid.sort_values(
        by=["cv_accuracy_mean", "simplicity", "n_comp_sort"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    return valid.iloc[0].to_dict(), table


def parse_window(value: Any) -> tuple[float, float]:
    if isinstance(value, str):
        value = literal_eval(value)
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise ValueError(f"Invalid window value: {value!r}.")
    return float(value[0]), float(value[1])


def evaluate_configuration(
    train_file: str | Path,
    test_file: str | Path,
    config: dict[str, Any],
) -> dict[str, Any]:
    method = str(config["method"])
    classifier = str(config["classifier"])
    window = parse_window(config["window"])
    n_components = config.get("n_components", FIXED_N_COMPONENTS)
    if pd.isna(n_components):
        n_components = FIXED_N_COMPONENTS
    n_components = int(n_components)

    train_blocks, y_train = load_epoch_blocks(train_file, method=method, window=window)
    test_blocks, y_test = load_epoch_blocks(test_file, method=method, window=window)
    x_train, weights = fit_features(
        train_blocks, y_train, method=method, n_components=n_components
    )
    x_test = transform_features(test_blocks, method=method, weights=weights)

    model = get_classifier(classifier)
    model.fit(x_train, y_train)
    predictions = model.predict(x_test)

    accuracy = float(accuracy_score(y_test, predictions))
    matrix = confusion_matrix(y_test, predictions, labels=LABEL_ORDER)
    n_correct = int(np.sum(predictions == y_test))
    n_test = int(y_test.size)
    p_value = float(binomtest(n_correct, n_test, p=0.5, alternative="greater").pvalue)

    return {
        "accuracy": accuracy,
        "accuracy_percent": accuracy * 100,
        "n_test_trials": n_test,
        "n_correct": n_correct,
        "p_vs_chance_0.5": p_value,
        "confusion_matrix": matrix.tolist(),
        "predictions": predictions,
        "y_test": y_test,
    }


def evaluate_fixed_csp_lda(train_file: str | Path, test_file: str | Path) -> dict[str, Any]:
    return evaluate_configuration(
        train_file,
        test_file,
        {
            "method": "csp",
            "classifier": "lda",
            "window": FIXED_WINDOW,
            "n_components": FIXED_N_COMPONENTS,
        },
    )
