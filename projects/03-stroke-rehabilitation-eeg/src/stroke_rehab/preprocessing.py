"""Continuous EEG preprocessing."""

from __future__ import annotations

import numpy as np
from scipy.signal import butter, filtfilt, iirnotch


def remove_dc(y: np.ndarray) -> np.ndarray:
    return y - np.mean(y, axis=0, keepdims=True)


def notch_filter(y: np.ndarray, fs: int, freq: float = 50.0, quality: float = 30.0) -> np.ndarray:
    nyquist = fs / 2.0
    if freq >= nyquist:
        return y
    b, a = iirnotch(w0=freq / nyquist, Q=quality)
    return filtfilt(b, a, y, axis=0)


def bandpass_filter(
    y: np.ndarray,
    fs: int,
    low: float = 8.0,
    high: float = 30.0,
    order: int = 4,
) -> np.ndarray:
    nyquist = fs / 2.0
    high = min(high, nyquist - 1e-6)
    if low <= 0 or high <= low:
        raise ValueError(f"Invalid band: low={low}, high={high}, fs={fs}.")
    b, a = butter(order, [low / nyquist, high / nyquist], btype="bandpass")
    return filtfilt(b, a, y, axis=0)


def preprocess_continuous(
    y: np.ndarray,
    fs: int,
    band: tuple[float, float] = (8.0, 30.0),
    use_notch: bool = True,
) -> np.ndarray:
    x = remove_dc(np.asarray(y, dtype=float))
    if use_notch:
        x = notch_filter(x, fs, freq=50.0, quality=30.0)
    return bandpass_filter(x, fs, low=band[0], high=band[1], order=4)
