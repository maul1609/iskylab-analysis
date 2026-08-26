"""Experiment timing/group metadata used by the iSKYLAB analysis scripts.

Times are seconds from the experiment time origin used in the supplied CSVs.
Keeping this information in one module avoids the previous duplication between
batch-running and mixing-analysis scripts.
"""

from __future__ import annotations

import numpy as np

CLOUD_WINDOWS = {
    "Exp006": np.array([14.0, 20.0]) * 60.0,
    "Exp007": np.array([18.0, 21.0]) * 60.0,
    "Exp008": np.array([15.0, 18.0]) * 60.0,
    "Exp009": np.array([13.0, 15.0]) * 60.0,
    "Exp010": np.array([13.0, 15.0]) * 60.0,
    "Exp011": np.array([17.0, 19.0]) * 60.0,
    "Exp012": np.array([14.0, 16.0]) * 60.0,
    "Exp013": np.array([13.0, 15.0]) * 60.0,
    "Exp014": np.array([15.0, 17.0]) * 60.0,
    "Exp015": np.array([10.5, 13.0]) * 60.0,
    "Exp016": np.array([10.5, 13.0]) * 60.0,
    "Exp017": np.array([1.67, 2.0]) * 60.0,
    "Exp018": np.array([9.0, 11.0]) * 60.0,
    "Exp019": np.array([14.0, 16.0]) * 60.0,
    "Exp028": np.array([11.0, 13.0]) * 60.0,
}

CLOUD_ONSET = {
    "Exp002": 5.0 * 60.0,
    "Exp003": 10.9 * 60.0,
    "Exp004": 5.8 * 60.0,
    "Exp005": 6.3 * 60.0,
    "Exp006": 9.0 * 60.0,
    "Exp007": 12.0 * 60.0,
    "Exp008": 10.0 * 60.0,
    "Exp009": 9.0 * 60.0,
    "Exp010": 9.0 * 60.0,
    "Exp011": 11.0 * 60.0,
    "Exp012": 10.0 * 60.0,
    "Exp013": 10.0 * 60.0,
    "Exp014": 9.0 * 60.0,
    "Exp015": 9.0 * 60.0,
    "Exp016": 9.0 * 60.0,
    "Exp017": 1.25 * 60.0,
    "Exp018": 7.0 * 60.0,
    "Exp019": 10.0 * 60.0,
    "Exp021": 12.0 * 60.0,
    "Exp022": 14.0 * 60.0,
    "Exp023": 13.0 * 60.0,
    "Exp025": 10.0 * 60.0,
    "Exp026": 12.0 * 60.0,
    "Exp027": 14.0 * 60.0,
    "Exp028": 7.0 * 60.0,
    "Exp029": 8.0 * 60.0,
}

# Historical sensitivity groups.  GROUP_TYPE determines which aerosol mode is
# shown on the x axis in the summary plot (1=first composition, 2=second).
BATCH_GROUPS = [
    ["Exp028", "Exp006", "Exp007"],
    ["Exp028", "Exp008", "Exp009"],
    ["Exp028", "Exp010"],
    ["Exp028", "Exp011"],
    ["Exp013", "Exp014"],
    ["Exp016", "Exp017"],
    ["Exp015", "Exp019"],
]
GROUP_TYPE = [2, 2, 2, 1, 1, 1, 2]
GROUP_UPDRAFT = [
    [1.3, 1.3, 1.3],
    [1.3, 1.3, 1.3],
    [1.3, 1.3],
    [1.3, 1.3],
    [1.3, 1.3],
    [1.3, 10.4],
    [1.3, 1.3],
]
