.. _plugins:

.. title:: ``z-brains`` plugins

Exploring ``z-brains`` plugins
============================================================

zbrains can accept data outside of micapipe and hippunfold as a 'plugin' folder. However, these data MUST be formatted as BIDS-derivatives exactly as in micapipe and hippunfold. If hippocampal surface data are present then they will be used but otherwise volumetric data will be mapped to hippocampal and subcortical surfaces. 

For example, the following plugin will import the data that is already mapped to hippocampal surfaces:

.. code-block:: console

    /data/mica3/BIDS_MICs/derivatives/plugin-InTS/
    └── sub-HC002
        └── ses-01
            ├── maps
            │   ├── sub-HC002_ses-01_hemi-L_space-T1w_den-2mm_label-hipp_InTS.shape.gii
            │   ├── sub-HC002_ses-01_hemi-L_surf-fsLR-32k_label-midthickness_InTS.func.gii
            │   ├── sub-HC002_ses-01_hemi-R_space-T1w_den-2mm_label-hipp_InTS.shape.gii
            │   ├── sub-HC002_ses-01_hemi-R_surf-fsLR-32k_label-midthickness_InTS.func.gii
            │   └── sub-HC002_ses-01_space-nativepro_map-InTS.nii.gz
            └── surf
                ├── sub-HC002_ses-01_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii
                ├── sub-HC002_ses-01_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii
                ├── sub-HC002_ses-01_hemi-L_space-T1w_den-2mm_label-hipp_midthickness.surf.gii
                └── sub-HC002_ses-01_hemi-R_space-T1w_den-2mm_label-hipp_midthickness.surf.gii

Similarly, the following input plugin will automatically map volumetric data to the corresponding surfaces:

.. code-block:: console

    /data/mica3/BIDS_MICs/derivatives/plugin-InTS/
    └── sub-HC002
        └── ses-01
            ├── maps
            │   ├── sub-HC002_ses-01_hemi-L_space-T1w_den-2mm_label-hipp_InTS.shape.gii
            │   ├── sub-HC002_ses-01_hemi-R_space-T1w_den-2mm_label-hipp_InTS.shape.gii
            │   └── sub-HC002_ses-01_space-nativepro_map-InTS.nii.gz
            └── surf
                ├── sub-HC002_ses-01_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii
                └── sub-HC002_ses-01_hemi-R_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii

For many calculations, its better to perform calculations on a surface rather than in a volume. Thus, when hippocampal surface data are available, they are prioritized. If they are not available, a volumetric image will be mapped to the surfaces in the `huippunfold` directory. 

In this example, we consdered the feature `InTS` or "intrinsic timescale". Thus we should add this to the list of features for zbrains to consider, but since it comes from a plugin, we add the prefix `plugin-`.

For example, zbrains can be run with the following command:

.. code-block:: console
    zbrains --sub HC002 --ses 01 --dataset /data/mica3/BIDS_MICs --zbrains OUTPUT --micapipe micapipe_v0.2.0  --hippunfold hippunfold_v1.3.0 --plugin plugin-INts --feat plugin-InTS

If we also want to consider the usual default zbrains features, we would run:

.. code-block:: console
    zbrains --sub HC002 --ses 01 --dataset /data/mica3/BIDS_MICs --zbrains OUTPUT --micapipe micapipe_v0.2.0  --hippunfold hippunfold_v1.3.0 --plugin plugin-INts --feat ADC FA flair qT1 thickness plugin-InTS

