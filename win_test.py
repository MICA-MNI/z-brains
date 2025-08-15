from zbrains.dataset import zbdataset, demographics
from zbrains.environment import zbenv
import os
import tempfile

features = ["FA", "ADC", "thickness", "qT1", "qT1-blur"]

env = zbenv(connectome_workbench_path="C:/Users/Ian/Downloads/workbench-windows64-v1.5.0/workbench/bin_windows64", num_threads=2, num_threads_wb=8)

control = demographics("data/participants_mics_hc_all.csv", normative_columns=["AGE", "SEX"], normative_dtypes=["int", "binary"])

control_dataset = zbdataset("controls", 
                    demographics=control,
                    micapipe_directory="E:/data/derivatives/micapipe", 
                    hippunfold_directory="E:/data/derivatives/hippunfold",
                    freesurfer_directory="E:/data/derivatives/freesurfer",
                    cortex=True,
                    hippocampus=True,
                    subcortical=True,
                    )

# control_dataset.process(output_directory="E:/zbrains_outputs", 
#                         features=features, 
#                         cortical_smoothing=10, 
#                         hippocampal_smoothing=5, 
#                         env=env, 
#                         verbose=True)

control_dataset.validate(output_directory="E:/zbrains_outputs", 
                         features=features, 
                         cortical_smoothing=10, 
                         hippocampal_smoothing=5, 
                         verbose=True)

patient = demographics("data/participants_mics_px_all.csv", reference=control, normative_columns=["AGE", "SEX"], normative_dtypes=["int", "binary"])
patient_dataset = zbdataset("patients", 
                    demographics=patient,
                    micapipe_directory="E:/data/derivatives/micapipe", 
                    hippunfold_directory="E:/data/derivatives/hippunfold",
                    freesurfer_directory="E:/data/derivatives/freesurfer",
                    cortex=True,
                    hippocampus=True,
                    subcortical=True,
                    )

# patient_dataset.process(output_directory="E:/zbrains_outputs", 
#                         features=features, 
#                         cortical_smoothing=10, 
#                         hippocampal_smoothing=5, 
#                         env=env, 
#                         verbose=True)


patient_dataset.validate(output_directory="E:/zbrains_outputs", 
                        features=features, 
                        cortical_smoothing=10, 
                        hippocampal_smoothing=5, 
                        verbose=True)


patient_dataset.analyze(output_directory="E:/zbrains_outputs", reference=control_dataset, method='wscore')

patient_dataset.clinical_report(
    output_directory="E:/zbrains_outputs", 
    approach='wscore',
    features=features,
    verbose=True
)

# patient_dataset.save_volumes(output_directory="/host/verges/tank/data/ian/zbrains_outputs", verbose=True)

# patient_dataset.save_DICOM(output_directory="/host/verges/tank/data/ian/zbrains_outputs", verbose=True)

print('end')