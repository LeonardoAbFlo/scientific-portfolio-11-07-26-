"""Log-variance, CSP and filter-bank CSP feature extraction."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.linalg import eigh

from .config import FBCSP_BANDS, FIXED_BAND, LABEL_ORDER
from .epochs import load_preprocess_epoch


def logvar_features(epochs: np.ndarray) -> np.ndarray:
    variance = np.var(epochs, axis=1)
    return np.log(variance + 1e-10)


def covariance_matrix(epoch: np.ndarray) -> np.ndarray:
    x = epoch.T
    covariance = x @ x.T
    trace = np.trace(covariance)
    if trace <= 0:
        trace = 1e-10
    return covariance / trace


def fit_csp(epochs: np.ndarray, labels: np.ndarray, n_components: int = 4) -> np.ndarray:
    labels = np.asarray(labels)
    class_covariances: list[np.ndarray] = []
    for class_value in LABEL_ORDER:
        class_epochs = epochs[labels == class_value]
        if class_epochs.size == 0:
            raise ValueError(f"No epochs found for class {class_value}.")
        covariance = np.mean([covariance_matrix(epoch) for epoch in class_epochs], axis=0)
        covariance += 1e-6 * np.eye(covariance.shape[0])
        class_covariances.append(covariance)

    eigenvalues, eigenvectors = eigh(
        class_covariances[0], class_covariances[0] + class_covariances[1]
    )
    order = np.argsort(eigenvalues)
    n_components = min(int(n_components), eigenvectors.shape[1])
    half = n_components // 2
    if n_components % 2 == 0:
        selected = np.r_[order[:half], order[-half:]]
    else:
        selected = np.r_[order[:half], order[-(half + 1):]]
    return eigenvectors[:, selected]


def transform_csp(epochs: np.ndarray, weights: np.ndarray) -> np.ndarray:
    features: list[np.ndarray] = []
    for epoch in epochs:
        projected = weights.T @ epoch.T
        variance = np.var(projected, axis=1)
        features.append(np.log(variance + 1e-10))
    return np.asarray(features)


def method_bands(method: str) -> tuple[tuple[float, float], ...]:
    if method in {"logvar", "csp"}:
        return (FIXED_BAND,)
    if method == "fbcsp":
        return tuple(FBCSP_BANDS)
    raise ValueError(f"Unknown feature method: {method}.")


def load_epoch_blocks(
    path: str | Path,
    method: str,
    window: tuple[float, float],
) -> tuple[list[np.ndarray], np.ndarray]:
    """Load one epoch array per analysis band and verify label alignment."""
    blocks: list[np.ndarray] = []
    reference_labels: np.ndarray | None = None
    for band in method_bands(method):
        _, epochs, labels = load_preprocess_epoch(path, window=window, band=band)
        if reference_labels is None:
            reference_labels = labels
        elif not np.array_equal(reference_labels, labels):
            raise ValueError("Trial labels changed across filter-bank bands.")
        blocks.append(epochs)
    if reference_labels is None:
        raise ValueError(f"No epoch blocks loaded from {path}.")
    return blocks, reference_labels


def fit_features(
    epoch_blocks: Sequence[np.ndarray],
    labels: np.ndarray,
    method: str,
    n_components: int,
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Fit data-driven transforms and return training features plus state."""
    if method == "logvar":
        return logvar_features(epoch_blocks[0]), []

    weights: list[np.ndarray] = []
    feature_blocks: list[np.ndarray] = []
    for epochs in epoch_blocks:
        csp_weights = fit_csp(epochs, labels, n_components=n_components)
        weights.append(csp_weights)
        feature_blocks.append(transform_csp(epochs, csp_weights))
    return np.concatenate(feature_blocks, axis=1), weights


def transform_features(
    epoch_blocks: Sequence[np.ndarray],
    method: str,
    weights: Sequence[np.ndarray],
) -> np.ndarray:
    """Apply fitted transforms to validation or test epochs."""
    if method == "logvar":
        return logvar_features(epoch_blocks[0])
    if len(epoch_blocks) != len(weights):
        raise ValueError("Number of epoch blocks and CSP weight matrices does not match.")
    return np.concatenate(
        [transform_csp(epochs, csp_weights) for epochs, csp_weights in zip(epoch_blocks, weights)],
        axis=1,
    )
