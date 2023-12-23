import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import nibabel as nib


LABELS = [26, 18, 11, 17, 13, 12, 10, 58, 54, 50, 53, 52, 51, 49]
STRUCTURES = ['Laccumb', 'Lamyg', 'Lcaud', 'Lhippo', 'Lpal', 'Lput', 'Lthal', 'Raccumb', 'Ramyg', 'Rcaud', 'Rhippo',
              'Rpal', 'Rput', 'Rthal']


parser = argparse.ArgumentParser(description='Subcortical mapping')

parser.add_argument('-id', '--subject_id', type=str, required=True, help='Subject ID')
parser.add_argument('-i', '--image', type=Path, default=None, help="Intensity image.")
parser.add_argument('-s', '--seg', type=Path, default=None, help="Subcortical segmentation.")
parser.add_argument('-v', '--vol', type=Path, default=None, help="Path to aseg.stats file.")
parser.add_argument('-o', '--output', type=Path, required=True, help="Output csv.")
parser.add_argument('--include_icv', action='store_true', default=False,
                    help='Include IntraCranial volume in csv.')

args = parser.parse_args()

if args.seg and args.image and args.vol is None:
    seg = np.asanyarray(nib.load(args.seg).dataobj)
    img = np.asanyarray(nib.load(args.image).dataobj)
    cols = STRUCTURES + ['ICV'] if args.include_icv else STRUCTURES
    df = pd.DataFrame(np.nan, columns=cols, index=pd.Index([args.subject_id], name='SubjID'))
    for i, k in enumerate(LABELS):
        df.loc[args.subject_id, STRUCTURES[i]] = img[seg == k].mean()

elif args.seg is None and args.image is None and args.vol:
    df_volumes = pd.read_csv(args.vol, comment="#", index_col=0, sep=r"\s+", header=None, usecols=[1, 3, 4],
                             names=['label', 'volume', 'structure'])
    volumes = df_volumes.loc[LABELS].volume.to_numpy()
    cols = STRUCTURES
    if args.include_icv:
        icv = float(next((line.split(',')[-2].strip() for line in open(args.vol)
                          if line.startswith("# Measure EstimatedTotalIntraCranialVol")), None))
        volumes = np.r_[volumes, icv]
        cols = STRUCTURES + ['ICV']

    df = pd.DataFrame(volumes[None], columns=cols, index=pd.Index([args.subject_id], name='SubjID'))

else:
    raise ValueError("Ambiguous inputs. Please provide either '--image' and '--seg' together or '--vol' alone.")

df.to_csv(args.output)
