.. _outputs:

.. title:: zbrains outputs

Exploring ``zbrains`` outputs
============================================================

What can you expect to find in your zbrains output directory? Let's investigate...


Raw features: maps
--------------------------------------------------------

Raw features sampled from neocortical and hippocampal surfaces, as well as subcortical structure volumes, can be found under the patients ``maps`` directory: ``<outputDirectory>/zbrains/<sub>/<ses>/maps``.
The ``maps`` directory thus contains three subdirectories:

.. parsed-literal::

    - <mapsdir>/cortex
    - <mapsdir>/hippocampus
    - <mapsdir>/subcortex

The outputs can be broadly categorized according to the region and modality they are sampling.

.. parsed-literal::

    - **Cortical outputs:**
        - <mapsdir>/cortex>/<BIDS_ID>_hemi-<X>_surf-fsLR-<X>k_feature-<X>_smooth-<X>mm.func.gii

    - **Hippocampal outputs:**
        - <mapsdir>/hippocampus>/<BIDS_ID>_hemi-<X>_den-<X>mm_feature-<X>_smooth-<X>mm.func.gii

    - **Subcortical outputs:**
        - <mapsdir>/subcortex>/<BIDS_ID>_feature-<X>.csv


Normalized features: Z-scored
--------------------------------------------------------

Z-scored features sampled from neocortical and hippocampal surfaces, as well as subcortical structure volumes, can be found under the patients ``norm-z`` directory: ``<outputDirectory>/zbrains/<sub>/<ses>/norm-z``.
The ``norm-z`` directory thus contains three subdirectories:

.. parsed-literal::

    - <norm-z>/cortex
    - <norm-z>/hippocampus
    - <norm-z>/subcortex

The outputs can be broadly categorized according to the region and modality they are sampling, as well as the analysis performed (regional vs. asymmetry) and the threshold defined by the user.

.. parsed-literal::

    - **Cortical outputs:**
        - <norm-z>/cortex>/<BIDS_ID>_hemi-<X>_surf-fsLR-<X>k_feature-<X>_smooth-<X>mm_<analysis>_<threshold>.func.gii

    - **Hippocampal outputs:**
        - <norm-z>/hippocampus>/<BIDS_ID>_hemi-<X>_den-<X>mm_feature-<X>_smooth-<X>mm_<analysis>_<threshold>.func.gii

    - **Subcortical outputs:**
        - <norm-z>/subcortex>/<BIDS_ID>_feature-<X>_<analysis>_<threshold>.csv


Normalized features: Normative model
--------------------------------------------------------

🚧 🚧 🚧 Under construction 🚧 🚧 🚧
