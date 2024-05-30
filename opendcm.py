import pydicom
import numpy as np

filepath = (
    "C:/Users/Ian/Documents/MR.1.3.12.2.1107.5.2.43.167017.20220920125221127268193"
)
fp = "E:/BIDS_MICS_Test/data/derivatives/Test/sub-PX001/ses-02/norm-z-volumetric/DICOM/sub-PX001_ses-02_surf-fsLR-32k_label-midthickness_feature-ADC_smooth-ctx-10mm_smooth-hipp-5mm_analysis-asymmetry/slice0000.dcm"
fp = "E:/BIDS_MICS_Test/data/derivatives/Test/sub-PX001/ses-02/norm-z-volumetric/DICOM/baseline/IM_0003.dcm"
fpreal = "E:/BIDS_MICS_Test/data/derivatives/Test/sub-PX001/ses-02/norm-z-volumetric/DICOM/baseline/IM_0001.dcm"
newfp = "C:/Users/Ian/Downloads/raw-RGB.dcm"


# ds = pydicom.filereader.dcmread(fp)
# dsreal = pydicom.filereader.dcmread(fpreal)
dsnew = pydicom.filereader.dcmread(newfp)
ar = dsnew.pixel_array
array = np.frombuffer(ar, dtype=np.uint8)
print("original ", dsnew)
# print("real ", dsreal)
