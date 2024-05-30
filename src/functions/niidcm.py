# Description: Convert NIfTI to DICOM
import os
from os.path import abspath
import numpy as np
import pydicom as pyd
from nii2dcm.dcm import DicomMRI

eps = 0.00000001192093


def write_slice(dcm, img_data, slice_index, output_dir):
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

    # img_slice = np.moveaxis(img_slice, -1, 0)
    # Instance UID – unique to current slice
    dcm.ds.SOPInstanceUID = pyd.uid.generate_uid(None)

    # write pixel data
    dcm.ds.PixelData = img_slice.tobytes()
    # write DICOM file
    dcm.ds.save_as(os.path.join(output_dir, output_filename), write_like_original=False)


def get_nii2dcm_parameters(nib_nii):
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

    # load nifti
    nii_img = nib_nii.get_fdata()

    # volume dimensions

    nX, nY, nZ, nF = (
        nib_nii.header["dim"][1],
        nib_nii.header["dim"][2],
        nib_nii.header["dim"][3],
        1,
    )
    dimX, dimY, dimZ = (
        nib_nii.header["pixdim"][1],
        nib_nii.header["pixdim"][2],
        nib_nii.header["pixdim"][3],
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
    A = nib_nii.affine
    dircosX = -1 * A[:3, 0] / dimX
    dircosY = -1 * A[:3, 1] / dimY

    image_pos_patient_array = []
    for iInstance in range(0, nInstances):
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
        # alternative:
        # 'ImageOrientationPatient': [dircosX[0], dircosX[1], dircosX[2], dircosY[0], dircosY[1], dircosY[2]],
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


def transfer_nii_hdr_series_tags(dcm, nii2dcm_parameters):
    """
    Transfer NIfTI header parameters applicable across Series

    dcm – nii2dcm DICOM object
    nii2dcm_parameters - parameters from NIfTI file
    """

    dcm.ds.Rows = nii2dcm_parameters["Rows"]
    dcm.ds.Columns = nii2dcm_parameters["Columns"]
    dcm.ds.PixelSpacing = [
        round(float(nii2dcm_parameters["dimX"]), 2),
        round(float(nii2dcm_parameters["dimY"]), 2),
    ]
    dcm.ds.SliceThickness = nii2dcm_parameters["SliceThickness"]
    dcm.ds.SpacingBetweenSlices = round(
        float(nii2dcm_parameters["SpacingBetweenSlices"]), 2
    )
    dcm.ds.ImageOrientationPatient = nii2dcm_parameters["ImageOrientationPatient"]
    dcm.ds.AcquisitionMatrix = nii2dcm_parameters["AcquisitionMatrix"]
    dcm.ds.SmallestImagePixelValue = (
        int(nii2dcm_parameters["SmallestImagePixelValue"])
        if int(nii2dcm_parameters["SmallestImagePixelValue"]) > 0
        else 0
    )  # SmallestImagePixelValue must be >= 0
    dcm.ds.LargestImagePixelValue = int(nii2dcm_parameters["LargestImagePixelValue"])
    dcm.ds.WindowCenter = nii2dcm_parameters["WindowCenter"]
    dcm.ds.WindowWidth = nii2dcm_parameters["WindowWidth"]
    dcm.ds.RescaleIntercept = nii2dcm_parameters["RescaleIntercept"]
    dcm.ds.RescaleSlope = nii2dcm_parameters["RescaleSlope"]
    dcm.file_meta.ImplementationVersionName = nii2dcm_parameters[
        "ImplementationVersionName"
    ]
    dcm.ds.AccessionNumber = nii2dcm_parameters["AccessionNumber"]
    dcm.ds.InstitutionName = nii2dcm_parameters["InstitutionName"]
    dcm.ds.PatientName = nii2dcm_parameters["PatientName"]


def transfer_nii_hdr_instance_tags(dcm, nii2dcm_parameters, instance_index):
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

    dcm.ds.InstanceNumber = nii2dcm_parameters["InstanceNumber"][instance_index]
    dcm.ds.SliceLocation = nii2dcm_parameters["SliceLocation"][instance_index]
    dcm.ds.ImagePositionPatient = [
        str(round(nii2dcm_parameters["ImagePositionPatient"][instance_index][0], 2)),
        str(round(nii2dcm_parameters["ImagePositionPatient"][instance_index][1], 2)),
        str(round(nii2dcm_parameters["ImagePositionPatient"][instance_index][2], 2)),
    ]
    dcm.ds.PhotometricInterpretation = "RGB"
    dcm.ds.SamplesPerPixel = 3
    dcm.ds.BitsAllocated = 8
    dcm.ds.BitsStored = 8
    dcm.ds.PlanarConfiguration = 1
    dcm.ds.HighBit = 7
    dcm.ds.add_new(0x00280006, "US", 0)


def convert_nifti_to_dicom(
    img, out_dir, feature, smooth_ctx, smooth_hipp, analysis, metadata=None
):
    array = img.get_fdata()

    # array = (array - np.min(array)) / (np.max(array) - np.min(array))
    # array = array * np.iinfo(np.int16).max
    # array = array.astype(np.int16)

    nii2dcm_parameters = get_nii2dcm_parameters(img)

    dicom = DicomMRI("nii2dcm_dicom_mri.dcm")
    if metadata is not None:
        pc = ""
        metadata = metadata.to_dict()
        print(metadata)
        if "participant_id" in metadata:
            dicom.ds.PatientID = str(metadata["participant_id"][0])
            del metadata["participant_id"]
        if "session_id" in metadata:
            dicom.ds.StudyID = str(metadata["session_id"][0])
            del metadata["session_id"]
        if "site" in metadata:
            dicom.ds.ManufacturerModelName = str(metadata["site"][0])
            del metadata["site"]
        if "age" in metadata:
            dicom.ds.PatientAge = "0" + str(int(metadata["age"][0])) + "Y"
            del metadata["age"]
        if "sex" in metadata:
            dicom.ds.PatientSex = str(metadata["sex"][0])
            del metadata["sex"]
        if "scan_date" in metadata:
            strs = str(metadata["scan_date"][0]).split(".")
            dicom.ds.StudyDate = strs[2] + strs[1] + strs[0]
            del metadata["scan_date"]
        if len(metadata) > 0:
            for key, value in metadata.items():
                pc += f"{key}: {value[0]}\n"
            dicom.ds.PatientComments = pc
    dicom.ds.StudyDescription = f"{feature.upper()}, {analysis}, Smoothing - Cort: {smooth_ctx}, Hipp: {smooth_hipp}"

    transfer_nii_hdr_series_tags(dicom, nii2dcm_parameters)
    for instance_index in range(nii2dcm_parameters["NumberOfInstances"]):

        # Transfer Instance tags
        transfer_nii_hdr_instance_tags(dicom, nii2dcm_parameters, instance_index)

        # Write slice
        write_slice(dicom, array, instance_index, out_dir)

    print(f"nii2dcm: DICOM files written to: {abspath(out_dir)}")  # TODO use logger
