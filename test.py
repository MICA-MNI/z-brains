from src.dataset import zbdataset, demographics
from src.environment import zbenv

features = ["FA", "ADC", "thickness", "qT1", "qT1-blur", "FLAIR", "FLAIR-blur"]

# features = ["FA", "ADC", "thickness", "qT1"]

env = zbenv(connectome_workbench_path="/usr/bin/", num_threads=12)

control = demographics("data/participants_mics_hc_all.csv", normative_columns=["AGE", "SEX"], normative_dtypes=["int", "binary"])

control_dataset = zbdataset("controls", 
                    demographics=control,
                    micapipe_directory="/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0", 
                    hippunfold_directory="/data/mica3/BIDS_MICs/derivatives/hippunfold_v1.3.0",
                    freesurfer_directory="/data/mica3/BIDS_MICs/derivatives/freesurfer",
                    cortex=True,
                    hippocampus=True,
                    subcortical=True,
                    )

control_dataset.process(output_directory="/host/verges/tank/data/ian/zbrains_outputs", 
                        features=features, 
                        cortical_smoothing=10, 
                        hippocampal_smoothing=5, 
                        env=env, 
                        verbose=True)

# control_dataset.validate(output_directory="/host/verges/tank/data/ian/zbrains_outputs", 
#                          features=["FA", "ADC", "thickness", "qT1", "qT1-blur"], 
#                          cortical_smoothing=5, 
#                          hippocampal_smoothing=2, 
#                          verbose=True)

patient = demographics("data/participants_mics_px_all.csv", reference=control, normative_columns=["AGE", "SEX"], normative_dtypes=["int", "binary"])
patient_dataset = zbdataset("patients", 
                    demographics=patient,
                    micapipe_directory="/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0", 
                    hippunfold_directory="/data/mica3/BIDS_MICs/derivatives/hippunfold_v1.3.0",
                    freesurfer_directory="/data/mica3/BIDS_MICs/derivatives/freesurfer",
                    cortex=True,
                    hippocampus=True,
                    subcortical=True,
                    )

patient_dataset.process(output_directory="/host/verges/tank/data/ian/zbrains_outputs", 
                        features=features, 
                        cortical_smoothing=10, 
                        hippocampal_smoothing=5, 
                        env=env, 
                        verbose=True)


# patient_dataset.validate(output_directory="/host/verges/tank/data/ian/zbrains_outputs", 
#                         features=["FA", "ADC", "thickness", "qT1", "qT1-blur"], 
#                         cortical_smoothing=5, 
#                         hippocampal_smoothing=2, 
#                         verbose=True)


patient_dataset.analyze(output_directory="/host/verges/tank/data/ian/zbrains_outputs", reference=control_dataset, method='wscore')

patient_dataset.clinical_report(
    output_directory="/host/verges/tank/data/ian/zbrains_outputs", 
    approach='wscore',
    analyses=['regional', 'assymmetry'],
    features=features,
    verbose=True
)

# patient_dataset.save_volumes(output_directory="/host/verges/tank/data/ian/zbrains_outputs", verbose=True)

# patient_dataset.save_DICOM(output_directory="/host/verges/tank/data/ian/zbrains_outputs", verbose=True)

print('end')