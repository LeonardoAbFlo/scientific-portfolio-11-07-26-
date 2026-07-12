from __future__ import annotations

import numpy as np

from stroke_rehab.epochs import extract_epochs_from_continuous, find_trigger_onsets
from stroke_rehab.features import fit_csp, logvar_features, transform_csp


def test_find_trigger_onsets_detects_blocks() -> None:
    trigger = np.array([0, 1, 1, 0, -1, -1, 0, 1, 0])
    onsets, labels = find_trigger_onsets(trigger)
    assert onsets.tolist() == [1, 4, 7]
    assert labels.tolist() == [1, -1, 1]


def test_epoch_extraction_shape() -> None:
    fs = 10
    eeg = np.random.default_rng(1).normal(size=(120, 4))
    trigger = np.zeros(120)
    trigger[[10, 60]] = [1, -1]
    epochs, labels = extract_epochs_from_continuous(eeg, trigger, fs, window=(0.0, 2.0))
    assert epochs.shape == (2, 20, 4)
    assert labels.tolist() == [1, -1]


def test_csp_feature_dimensions() -> None:
    rng = np.random.default_rng(2)
    epochs = rng.normal(size=(20, 100, 8))
    labels = np.array([1] * 10 + [-1] * 10)
    weights = fit_csp(epochs, labels, n_components=4)
    features = transform_csp(epochs, weights)
    assert weights.shape == (8, 4)
    assert features.shape == (20, 4)
    assert logvar_features(epochs).shape == (20, 8)
