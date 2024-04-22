import numpy as np
import pandas as pd
import nibabel as nib
from pathlib import Path



def subcortical_mapping(subject_id: str, image: Path = None, seg: Path = None, vol: Path = None, output: Path = None, include_icv: bool = False):
    LABELS = [26, 18, 11, 17, 13, 12, 10, 58, 54, 50, 53, 52, 51, 49]
    STRUCTURES = ['Laccumb', 'Lamyg', 'Lcaud', 'Lhippo', 'Lpal', 'Lput', 'Lthal', 'Raccumb', 'Ramyg', 'Rcaud', 'Rhippo',
                'Rpal', 'Rput', 'Rthal']
    if seg and image and vol is None:
        seg = np.asanyarray(nib.load(seg).dataobj)
        img = np.asanyarray(nib.load(image).dataobj)
        cols = STRUCTURES + ['ICV'] if include_icv else STRUCTURES
        df = pd.DataFrame(np.nan, columns=cols, index=pd.Index([subject_id], name='SubjID'))
        for i, k in enumerate(LABELS):
            df.loc[subject_id, STRUCTURES[i]] = img[seg == k].mean()

    elif seg is None and image is None and vol:
        df_volumes = pd.read_csv(vol, comment="#", index_col=0, sep=r"\s+", header=None, usecols=[1, 3, 4],
                                 names=['label', 'volume', 'structure'])
        volumes = df_volumes.loc[LABELS].volume.to_numpy()
        cols = STRUCTURES
        if include_icv:
            icv = float(next((line.split(',')[-2].strip() for line in open(vol)
                              if line.startswith("# Measure EstimatedTotalIntraCranialVol")), None))
            volumes = np.r_[volumes, icv]
            cols = STRUCTURES + ['ICV']

        df = pd.DataFrame(volumes[None], columns=cols, index=pd.Index([subject_id], name='SubjID'))

    else:
        raise ValueError("Ambiguous inputs. Please provide either 'image' and 'seg' together or 'vol' alone.")

    df.to_csv(output)