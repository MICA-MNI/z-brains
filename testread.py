from dataclasses import dataclass
from operator import itemgetter
import numpy as np
import nibabel as nb
from typing import Tuple


@dataclass
class VolGeom:
    # Using tuples to show shape. Will mostly be numpy arrays.
    shape: Tuple[int, int, int]
    zooms: Tuple[float, float, float]
    cosines: Tuple[
        Tuple[float, float, float],
        Tuple[float, float, float],
        Tuple[float, float, float],
    ]
    c_ras: Tuple[float, float, float]

    @classmethod
    def from_mghheader(cls, header):
        return cls(
            shape=header.dims[:3],
            zooms=header.delta,
            cosines=header.Mdc.T,
            c_ras=header.Pxyz_c,
        )

    @classmethod
    def from_surf_footer(cls, footer):
        return cls(
            shape=footer["volume"],
            zooms=footer["voxelsize"],
            cosines=np.vstack(itemgetter("xras", "yras", "zras")(footer)).T,
            c_ras=footer["cras"],
        )

    @classmethod
    def from_gifti_metadata(cls, meta):
        return cls(
            shape=np.int16(
                (meta["VolGeomWidth"], meta["VolGeomHeight"], meta["VolGeomDepth"])
            ),
            zooms=np.float64(
                (meta["VolGeomXsize"], meta["VolGeomYsize"], meta["VolGeomZsize"])
            ),
            cosines=np.float64(
                itemgetter(*(f"VolGeom{col}_{row}" for row in "RAS" for col in "XYZ"))(
                    meta
                )
            ).reshape((3, 3)),
            c_ras=np.float64(
                (meta["VolGeomC_R"], meta["VolGeomC_A"], meta["VolGeomC_S"])
            ),
        )

    def tkreg_affine(self):
        tkrcosines = np.array([[-1, 0, 0], [0, 0, 1], [0, -1, 0]])
        mat = tkrcosines * self.zooms
        return nb.affines.from_matvec(mat, -mat @ self.shape / 2)

    def scanner_affine(self):
        mat = self.cosines * self.zooms
        return nb.affines.from_matvec(mat, self.c_ras - mat @ self.shape / 2)

    def tkreg2scanner(self):
        return self.scanner_affine() @ np.linalg.inv(self.tkreg_affine())


orig = nb.load(
    "E:/BIDS_MICS_Test/data/derivatives/micapipe/sub-HC005/ses-01/surf/sub-HC005_ses-01_hemi-L_surf-fsnative_label-midthickness.surf.gii"
)
new = nb.load(
    "E:/data/derivatives/micapipe/sub-PX001/ses-02/surf/sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii"
)
coords, triangles = orig.darrays
geom = VolGeom.from_gifti_metadata(coords.meta)

coords, triangles = new.darrays
scanner_coords = coords.__class__(
    nb.affines.apply_affine(geom.tkreg2scanner(), coords.data).astype("f4"),
    meta=coords.meta,
)
ourscanner = nb.GiftiImage(darrays=[scanner_coords, triangles], meta=new.meta)
ourscanner.to_filename(
    "C:/Users/Ian/Documents/GitHub/z-brains-IanTesting/src/data/input-white.surf.gii",
)
print(ourscanner)
