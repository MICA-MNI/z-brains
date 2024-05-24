import SimpleITK as sitk
import os
import time
import numpy as np


def convert_nifti_to_dicom(in_dir, out_dir):
    new_img = sitk.ReadImage(in_dir)

    array = sitk.GetArrayFromImage(new_img)
    # Normalize the array to the range [0, 1]
    array = (array - np.min(array)) / (np.max(array) - np.min(array))

    # Scale the array to the range [0, sitk.sitkInt16]
    array = array * np.iinfo(np.int16).max

    # Convert the array back to a SimpleITK image
    new_img = sitk.GetImageFromArray(array.astype(np.int16))
    print(array)
    modification_time = time.strftime("%H%M%S")
    modification_date = time.strftime("%Y%m%d")

    direction = new_img.GetDirection()
    series_tag_values = [
        ("0008|0031", modification_time),  # Series Time
        ("0008|0021", modification_date),  # Series Date
        ("0008|0008", "DERIVED\\SECONDARY"),  # Image Type
        (
            "0020|000e",
            "1.2.826.0.1.3680043.2.1125."
            + modification_date
            + ".1"
            + modification_time,
        ),  # Series Instance UID
        (
            "0020|0037",
            "\\".join(
                map(
                    str,
                    (
                        direction[0],
                        direction[3],
                        direction[6],  # Image Orientation (Patient)
                        direction[1],
                        direction[4],
                        direction[7],
                    ),
                )
            ),
        ),
        ("0008|103e", "Created-SimpleITK"),
    ]  # Series Description

    # Write slices to output directory
    list(
        map(
            lambda i: writeSlices(series_tag_values, new_img, i, out_dir),
            range(new_img.GetDepth()),
        )
    )


def writeSlices(series_tag_values, new_img, i, out_dir):

    image_slice = new_img[:, :, i]
    writer = sitk.ImageFileWriter()

    writer.KeepOriginalImageUIDOn()

    # Tags shared by the series.
    list(
        map(
            lambda tag_value: image_slice.SetMetaData(tag_value[0], tag_value[1]),
            series_tag_values,
        )
    )

    # Slice specific tags.
    image_slice.SetMetaData(
        "0008|0012", time.strftime("%Y%m%d")
    )  # Instance Creation Date
    image_slice.SetMetaData(
        "0008|0013", time.strftime("%H%M%S")
    )  # Instance Creation Time

    # Setting the type to CT preserves the slice location.
    image_slice.SetMetaData(
        "0008|0060", "CT"
    )  # set the type to CT so the thickness is carried over

    # (0020, 0032) image position patient determines the 3D spacing between slices.
    image_slice.SetMetaData(
        "0020|0032",
        "\\".join(map(str, new_img.TransformIndexToPhysicalPoint((0, 0, i)))),
    )  # Image Position (Patient)
    image_slice.SetMetaData("0020,0013", str(i))  # Instance Number

    # Write to the output directory and add the extension dcm, to force writing in DICOM format.
    writer.SetFileName(os.path.join(out_dir, f"slice{str(i).zfill(4)}.dcm"))
    writer.Execute(image_slice)
