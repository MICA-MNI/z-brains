import pandas as pd
import numpy as np
import nibabel as nib

SUBCORTICAL_LABELS = [26, 18, 11, 17, 13, 12, 10, 58, 54, 50, 53, 52, 51, 49]
SUBCORTICAL_STRUCTURES = [
    "Laccumb", "Lamyg", "Lcaud", "Lhippo", "Lpal", "Lput", "Lthal",
    "Raccumb", "Ramyg", "Rcaud", "Rhippo", "Rpal", "Rput", "Rthal",
]

# This script is now intended to be run as a Snakemake script, using snakemake.input and snakemake.output
def main():
    seg_file = snakemake.input["seg_file"]
    feature_map = snakemake.input["feature_map"]
    output_csv = snakemake.output["out"]

    seg = np.asarray(nib.load(seg_file).dataobj)
    img = np.asarray(nib.load(feature_map).dataobj)
    df = pd.DataFrame(np.nan, columns=SUBCORTICAL_STRUCTURES, index=pd.Index(["subj"], name="SubjID"))
    for i, label in enumerate(SUBCORTICAL_LABELS):
        region_mask = seg == label
        if np.any(region_mask):
            df.loc["subj", SUBCORTICAL_STRUCTURES[i]] = np.nanmean(img[region_mask])
    df.to_csv(output_csv)

if __name__ == "__main__":
    main() 