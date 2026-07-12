"""Central analysis settings inherited from the original notebook."""

from __future__ import annotations

SUBJECTS = ("P1", "P2", "P3")
STAGES = ("pre", "post")
RUN_TYPES = ("training", "test")

LABEL_ORDER = (1, -1)
LABEL_NAMES = ("Left imagery (+1)", "Right imagery (-1)")

FIXED_WINDOW = (2.0, 8.0)
FIXED_BAND = (8.0, 30.0)
FIXED_N_COMPONENTS = 4

WINDOW_GRID = (
    (2.0, 8.0),
    (2.0, 7.0),
    (3.0, 7.0),
    (4.0, 8.0),
)
METHOD_GRID = ("logvar", "csp", "fbcsp")
CLASSIFIER_GRID = ("lda", "svm")
N_COMPONENTS_GRID = (2, 4, 6)

FBCSP_BANDS = (
    (8.0, 12.0),
    (12.0, 16.0),
    (16.0, 20.0),
    (20.0, 24.0),
    (24.0, 30.0),
)

# Channel order used in the notebook's activity/time-course summaries.
ACTIVITY_CHANNELS = (
    "FC5", "FC1", "FCz", "FC2", "FC6",
    "C5", "C3", "C1", "Cz", "C2", "C4", "C6",
    "CP5", "CP1", "CP2", "CP6",
)

# Channel order used in the notebook's permutation-contribution analysis.
# Confirm this order against the acquisition metadata before interpretation.
CONTRIBUTION_CHANNELS = (
    "FC3", "FCz", "FC4",
    "C5", "C3", "C1", "Cz", "C2", "C4", "C6",
    "CP3", "CP1", "CPz", "CP2", "CP4",
    "Pz",
)

RANDOM_STATE = 42
