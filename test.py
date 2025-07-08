from src.dataset import zbdataset, demographics
from src.environment import zbenv



control = demographics("data/participants_mics_hc_all.csv", normative_columns=["AGE", "SEX"])

control_dataset = zbdataset("controls", 
                    demographics=control,
                    micapipe_directory="/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0", 
                    hippunfold_directory="/data/mica3/BIDS_MICs/derivatives/hippunfold_v1.3.0",
                    freesurfer_directory="/data/mica3/BIDS_MICs/derivatives/freesurfer",
                    cortex=True,
                    hippocampus=True,
                    subcortical=True,
                    )

control_dataset.add_features("FA", "ADC", "thickness", "FLAIR", "qT1", "FLAIR-blur", "qT1-blur", verbose=False)

env = zbenv(connectome_workbench_path="/usr/bin/")

control_dataset.process(output_directory="/host/verges/tank/data/ian/zbrains_outputs", cortical_smoothing=5, hippocampal_smoothing=2, env=env, verbose=True)

print('end')