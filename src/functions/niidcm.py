# Description: Convert NIfTI to DICOM
import os
from os.path import abspath
import numpy as np
import pydicom as pyd
from .dcmclasses import DicomMRI
from multiprocessing import Pool

eps = 0.00000001192093


def write_slice_wrapper(args):
    ds, array, instance_index, out_dir, nii2dcm_parameters = args
    transfer_nii_hdr_instance_tags(ds, nii2dcm_parameters, instance_index)
    write_slice(ds, array, instance_index, out_dir)


def write_slice(ds, img_data, slice_index, output_dir):
    """
    write a single DICOM slice

    dcm – nii2dcm DICOM object
    img_data - [nX, nY, nSlice] image pixel data, such as from NIfTI file
    slice_index – slice index in nibabel img_data array (important: counts from 0, whereas DICOM instances count from 1)
    output_dir – output DICOM file save location
    """

    output_filename = r"IM_%04d.dcm" % (
        slice_index + 1
    )  # begin filename from 1, e.g. IM_0001.dcm

    img_slice = img_data[:, :, slice_index]
    img_slice = img_slice.flatten().astype(np.uint8)

    # Instance UID – unique to current slice
    ds.SOPInstanceUID = pyd.uid.generate_uid(None)

    # write pixel data
    ds.PixelData = img_slice.tobytes()
    # write DICOM file
    ds.save_as(os.path.join(output_dir, output_filename), write_like_original=False)


def get_nii2dcm_parameters(header, affine, nii_img):
    """
    Get general NIfTI header parameters relevant for DICOM tag transferal.
    :nib_nii - NIfTI loaded with nibabel
    :nii_parameters - parameters to transfer to DICOM header
    """

    def fnT1N(A, N):
        # Subfn: calculate T1N vector
        # A = affine matrix [4x4]
        # N = slice number (counting from 1)
        T1N = A.dot([0, 0, N - 1, 1])
        return T1N

    # volume dimensions

    nX, nY, nZ, nF = (
        header["dim"][1],
        header["dim"][2],
        header["dim"][3],
        1,
    )
    dimX, dimY, dimZ = (
        header["pixdim"][1],
        header["pixdim"][2],
        header["pixdim"][3],
    )

    # Instances & Slice Spacing
    nInstances = nZ * nF
    sliceIndices = np.repeat(range(1, nZ + 1), nF)
    voxelSpacing = dimZ
    zLocLast = (voxelSpacing * nZ) - voxelSpacing
    sliceLoca = np.repeat(np.linspace(0, zLocLast, num=nZ), nF)

    # Windowing & Signal Intensity
    maxI = np.amax(nii_img)
    minI = np.amin(nii_img)
    windowCenter = round((maxI - minI) / 2)
    windowWidth = round(maxI - minI)
    rescaleIntercept = 0
    rescaleSlope = 1

    # FOV
    fovX = nX * dimX
    fovY = nY * dimY
    fovZ = nZ * dimZ

    # slice positioning in 3-D space
    # nb: -1 for dir cosines gives consistent orientation between Nifti and DICOM in ITK-Snap
    A = affine
    dircosX = -1 * A[:3, 0] / dimX
    dircosY = -1 * A[:3, 1] / dimY

    image_pos_patient_array = []
    for iInstance in range(nInstances):
        T1N = fnT1N(A, iInstance)
        image_pos_patient_array.append([T1N[0], T1N[1], T1N[2]])

    # output dictionary
    nii2dcm_parameters = {
        # series parameters
        "dimX": dimX,
        "dimY": dimY,
        "SliceThickness": str(dimZ),
        "SpacingBetweenSlices": str(dimZ),
        "AcquisitionMatrix": [0, nX, nY, 0],
        "Rows": nX,
        "Columns": nY,
        "NumberOfSlices": nZ,
        "NumberOfInstances": nZ * nF,
        "PixelSpacing": [dimX, dimY],
        "FOV": [fovX, fovY, fovZ],
        "SmallestImagePixelValue": minI,
        "LargestImagePixelValue": maxI,
        "WindowCenter": str(windowCenter),
        "WindowWidth": str(windowWidth),
        "RescaleIntercept": str(rescaleIntercept),
        "RescaleSlope": str(rescaleSlope),
        "SpacingBetweenSlices": round(float(dimZ), 2),
        "ImageOrientationPatient": [
            dircosY[0],
            dircosY[1],
            dircosY[2],
            dircosX[0],
            dircosX[1],
            dircosX[2],
        ],
        # instance parameters
        "InstanceNumber": sliceIndices,
        "SliceLocation": sliceLoca,
        "ImagePositionPatient": image_pos_patient_array,
        "ImplementationVersionName": "z-Brains v0.1.0",
        "AccessionNumber": "",
        "InstitutionName": "z-Brains Software",
        "PatientName": "",
    }

    return nii2dcm_parameters


def transfer_nii_hdr_series_tags(ds, nii2dcm_parameters, file_meta):
    """
    Transfer NIfTI header parameters applicable across Series

    dcm – nii2dcm DICOM object
    nii2dcm_parameters - parameters from NIfTI file
    """

    ds.Rows = nii2dcm_parameters["Rows"]
    ds.Columns = nii2dcm_parameters["Columns"]
    ds.PixelSpacing = [
        round(float(nii2dcm_parameters["dimX"]), 2),
        round(float(nii2dcm_parameters["dimY"]), 2),
    ]
    ds.SliceThickness = nii2dcm_parameters["SliceThickness"]
    ds.SpacingBetweenSlices = round(
        float(nii2dcm_parameters["SpacingBetweenSlices"]), 2
    )
    ds.ImageOrientationPatient = nii2dcm_parameters["ImageOrientationPatient"]
    ds.AcquisitionMatrix = nii2dcm_parameters["AcquisitionMatrix"]
    ds.SmallestImagePixelValue = max(
        int(nii2dcm_parameters["SmallestImagePixelValue"]), 0
    )
    ds.LargestImagePixelValue = int(nii2dcm_parameters["LargestImagePixelValue"])
    ds.WindowCenter = nii2dcm_parameters["WindowCenter"]
    ds.WindowWidth = nii2dcm_parameters["WindowWidth"]
    ds.RescaleIntercept = nii2dcm_parameters["RescaleIntercept"]
    ds.RescaleSlope = nii2dcm_parameters["RescaleSlope"]
    file_meta.ImplementationVersionName = nii2dcm_parameters[
        "ImplementationVersionName"
    ]
    ds.AccessionNumber = nii2dcm_parameters["AccessionNumber"]
    ds.InstitutionName = nii2dcm_parameters["InstitutionName"]
    ds.PatientName = nii2dcm_parameters["PatientName"]


def transfer_nii_hdr_instance_tags(ds, nii2dcm_parameters, instance_index):
    """
    Transfer NIfTI header parameters applicable to Instance

    dcm – nii2dcm DICOM object
    nii2dcm_parameters - parameters from NIfTI file
    instance_index - slice number in NIfTI file
    """

    # Possible per Instance Tags
    # SOPInstanceUID (set within write_slice function)
    # InstanceNumber
    # SliceLocation
    # ImagePositionPatient

    ds.InstanceNumber = nii2dcm_parameters["InstanceNumber"][instance_index]
    ds.SliceLocation = nii2dcm_parameters["SliceLocation"][instance_index]
    ds.ImagePositionPatient = [
        str(round(nii2dcm_parameters["ImagePositionPatient"][instance_index][0], 2)),
        str(round(nii2dcm_parameters["ImagePositionPatient"][instance_index][1], 2)),
        str(round(nii2dcm_parameters["ImagePositionPatient"][instance_index][2], 2)),
    ]
    ds.PhotometricInterpretation = "RGB"
    ds.SamplesPerPixel = 3
    ds.BitsAllocated = 8
    ds.BitsStored = 8
    ds.PlanarConfiguration = 1
    ds.HighBit = 7
    ds.add_new(0x00280006, "US", 0)


def convert_nifti_to_dicom(
    array,
    header,
    affine,
    out_dir,
    feature,
    smooth_ctx,
    smooth_hipp,
    analysis,
    metadata=None,
):
    """
    Convert a NIfTI image to DICOM format with specified metadata.

    This function converts a NIfTI image to DICOM format, incorporating metadata information
    such as participant details, study ID, age, sex, scan date, and comments.

    Args:
        array: NIfTI image array to convert to DICOM.
        out_dir: Output directory to save the DICOM files.
        feature: Specific feature of the image.
        smooth_ctx: Level of smoothing for the cortex.
        smooth_hipp: Level of smoothing for the hippocampus.
        analysis: Type of analysis performed on the image.
        metadata: Metadata information for the DICOM file (default is None).

    Returns:
        None
    """
    print("Starting conversion of NIfTI image to DICOM format...")
    nii2dcm_parameters = get_nii2dcm_parameters(header, affine, array)
    print("NIfTI to DICOM parameters obtained.")

    dicom = DicomMRI("nii2dcm_dicom_mri.dcm")
    ds = dicom.ds
    filemeta = dicom.file_meta
    print("DICOM object initialized.")
    if metadata is not None:
        pc = ""
        metadata = metadata.to_dict()
        if "participant_id" in metadata:
            ds.PatientID = str(metadata["participant_id"][0])
            del metadata["participant_id"]
        if "session_id" in metadata:
            ds.StudyID = str(metadata["session_id"][0])
            del metadata["session_id"]
        if "site" in metadata:
            ds.ManufacturerModelName = str(metadata["site"][0])
            del metadata["site"]
        if "age" in metadata:
            ds.PatientAge = "0" + str(int(metadata["age"][0])) + "Y"
            del metadata["age"]
        if "sex" in metadata:
            ds.PatientSex = str(metadata["sex"][0])
            del metadata["sex"]
        if "scan_date" in metadata:
            strs = str(metadata["scan_date"][0]).split(".")
            ds.StudyDate = strs[2] + strs[1] + strs[0]
            del metadata["scan_date"]
        if len(metadata) > 0:
            for key, value in metadata.items():
                pc += f"{key}: {value[0]}\n"
            ds.PatientComments = pc
    print("Metadata processed and transferred to DICOM object.")
    ds.StudyDescription = (
        f"{feature.upper()}, {analysis}, Smooth-Cort: {smooth_ctx}, Hipp: {smooth_hipp}"
    )

    print("Study description set.")

    transfer_nii_hdr_series_tags(ds, nii2dcm_parameters, filemeta)
    print("NIfTI header series tags transferred.")

    with Pool() as p:
        p.map(
            write_slice_wrapper,
            [
                (ds, array, i, out_dir, nii2dcm_parameters)
                for i in range(nii2dcm_parameters["NumberOfInstances"])
            ],
        )
    print("DICOM slices written.")

    print(f"nii2dcm: DICOM files written to: {abspath(out_dir)}")
    print("Conversion complete.")
