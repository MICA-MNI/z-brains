from typing import List, Literal

VERSION = "v0.1.0 'reborn-again-revengeance'"

Feature = Literal[
    "ADC", "FA", "flair", "qT1", "thickness", "qT1_blur", "flair_blur", "plugin_*"
]
LIST_FEATURES: List[Feature] = list(Feature.__args__)

map_feature_to_file = dict(
    thickness="thickness",
    volume="volume",
    flair="flair",
    ADC="ADC",
    FA="FA",
    qT1="T1map",
    qT1_blur="T1map_blur",
    flair_blur="flair_blur",
)

Structure = Literal["cortex", "subcortex", "hippocampus"]
LIST_STRUCTURES: List[str] = list(Structure.__args__)

Task = Literal["proc", "analysis"]
LIST_TASKS: List[str] = list(Task.__args__)

Approach = Literal["zscore", "norm"]
LIST_APPROACHES: List[str] = list(Approach.__args__)

Analysis = Literal["regional", "asymmetry"]
LIST_ANALYSES: List[str] = list(Analysis.__args__)

Resolution = Literal["low", "high"]
# Resolution = Literal['high']
LIST_RESOLUTIONS: List[str] = list(Resolution.__args__)

LIST_LABELS_CTX = ["white", "midthickness", "pial", "swm[0-9]+"]
LIST_LABELS_HIP = ["midthickness"]

# Subject dir folder names
FOLDER_LOGS = "logs"
FOLDER_MAPS = "maps"
FOLDER_NORM_Z = "norm-z"
FOLDER_NORM_MODEL = "norm-normative"

approach_to_folder = dict(zip(LIST_APPROACHES, [FOLDER_NORM_Z, FOLDER_NORM_MODEL]))

# Structure folders
FOLDER_CTX = "cortex"
FOLDER_SCTX = "subcortex"
FOLDER_HIP = "hippocampus"

struct_to_folder = dict(zip(LIST_STRUCTURES, [FOLDER_CTX, FOLDER_SCTX, FOLDER_HIP]))

# Default smoothing - in mm
DEFAULT_SMOOTH_CTX = 5
DEFAULT_SMOOTH_HIP = 2

# Default threshold
DEFAULT_THRESHOLD = 1.96

# Resolution
LOW_RESOLUTION_CTX = "5k"
HIGH_RESOLUTION_CTX = "32k"
LOW_RESOLUTION_HIP = "0p5mm"
HIGH_RESOLUTION_HIP = "2mm"

map_resolution_ctx = {"low": LOW_RESOLUTION_CTX, "high": HIGH_RESOLUTION_CTX}
map_resolution_hip = {"low": LOW_RESOLUTION_HIP, "high": HIGH_RESOLUTION_HIP}


class ProcessingException(Exception):
    pass
