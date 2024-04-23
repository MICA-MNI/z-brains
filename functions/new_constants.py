LIST_TASKS = ['proc', 'analysis']
LIST_STRUCTURES = ['cortex', 'subcortex', 'hippocampus']
LIST_FEATURES = ['ADC', 'FA', 'flair', 'T1map', 'thickness']
LIST_RESOLUTIONS = ['low', 'high']
DEFAULT_SMOOTH_CTX = 5
DEFAULT_SMOOTH_HIP = 2
DEFAULT_THRESHOLD = 1.96
FOLDER_LOGS = 'logs'
FOLDER_MAPS = 'maps'
FOLDER_NORM_Z = 'norm-z'
FOLDER_NORM_MODEL = 'norm-normative'
FOLDER_CTX = 'cortex'
FOLDER_SCTX = 'subcortex'
FOLDER_HIP = 'hippocampus'
LOW_RESOLUTION_CTX = '5k'
HIGH_RESOLUTION_CTX = '32k'
LOW_RESOLUTION_HIP = '0p5mm'
HIGH_RESOLUTION_HIP = '2mm'
LIST_APPROACHES = ['zscore', 'norm']
LIST_ANALYSES = ['regional', 'asymmetry']
Resolution = ['low', 'high']
Resolution = ['high']
map_resolution_ctx = {'low': LOW_RESOLUTION_CTX, 'high': HIGH_RESOLUTION_CTX}
map_resolution_hip = {'low': LOW_RESOLUTION_HIP, 'high': HIGH_RESOLUTION_HIP}

class ProcessingException(Exception):
    pass
