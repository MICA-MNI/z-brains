
#!/usr/bin/env python3
import pickle
import numpy as np
import pandas as pd
import nibabel as nib
from pathlib import Path

def load_gifti(path):
    img = nib.load(path)
    arrays = [da.data for da in img.darrays]
    arr = np.stack(arrays, axis=0) if len(arrays) > 1 else arrays[0]
    # if it’s >2D, collapse extras
    if arr.ndim > 2:
        arr = arr.mean(axis=tuple(range(2, arr.ndim)))
    return arr

data = np.stack([load_gifti(str(p)) for p in snakemake.input.giis], axis=0)

# save as a pickle
with open(snakemake.output.mat, "wb") as f:
    pickle.dump(data, f)
