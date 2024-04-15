from typing import List

VERSION = "v0.0.2 'reborn'"

Feature = ['ADC', 'FA', 'flair', 'qT1', 'thickness']
LIST_FEATURES: List[str] = Feature

map_feature_to_file = dict(thickness='thickness', volume='volume',
                           flair='flair', ADC='ADC', FA='FA', qT1='T1map')

Structure = ['cortex', 'subcortex', 'hippocampus']
LIST_STRUCTURES: List[str] = Structure

Task = ['proc', 'analysis']
LIST_TASKS: List[str] = Task

Approach = ['zscore', 'norm']
LIST_APPROACHES: List[str] = Approach

Analysis = ['regional', 'asymmetry']
LIST_ANALYSES: List[str] = Analysis

Resolution = ['low', 'high']
Resolution = ['high']
LIST_RESOLUTIONS: List[str] = Resolution

LIST_LABELS_CTX = ['white', 'midthickness', 'pial', 'swm[0-9]+']
LIST_LABELS_HIP = ['midthickness']

# Subject dir folder names
FOLDER_LOGS = 'logs'
FOLDER_MAPS = 'maps'
FOLDER_NORM_Z = 'norm-z'
FOLDER_NORM_MODEL = 'norm-normative'

approach_to_folder = dict(zip(LIST_APPROACHES,
                              [FOLDER_NORM_Z, FOLDER_NORM_MODEL]))

# Structure folders
FOLDER_CTX = 'cortex'
FOLDER_SCTX = 'subcortex'
FOLDER_HIP = 'hippocampus'

struct_to_folder = dict(zip(LIST_STRUCTURES,
                            [FOLDER_CTX, FOLDER_SCTX, FOLDER_HIP]))

# Default smoothing - in mm
DEFAULT_SMOOTH_CTX = 5
DEFAULT_SMOOTH_HIP = 2

# Default threshold
DEFAULT_THRESHOLD = 1.96

# Resolution
LOW_RESOLUTION_CTX = '5k'
HIGH_RESOLUTION_CTX = '32k'
LOW_RESOLUTION_HIP = '0p5mm'
HIGH_RESOLUTION_HIP = '2mm'


# Do not change - this is used in bash
if __name__ == '__main__':
    constants = {k: v for k, v in globals().items() if k.isupper()}
    for k, v in constants.items():
        if isinstance(v, list):
            print('export {}="{}"'.format(k, " ".join(v)))
        else:
            print('export {}="{}"'.format(k, v))