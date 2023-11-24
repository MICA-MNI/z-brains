.. _tutinstruc:

.. title:: Getting started

Usage notes
======================================================================

Supported imaging modalities
--------------------------------------------------------

The ``z-brains`` software  currently supports four MRI modalities and associated feature maps for normalization:

.. list-table::
  :widths: 10 1000
  :header-rows: 1

  * - **Modalities**
    - **Features**
  * - T1-weighted (T1w)
    - Vertexwise neocortical & hippocampal thickness, and average volume in 14 subcortical structures
  * - Quantitative T1 (qT1)
    - Vertexwise neocortical & hippocampal T1 relaxation time, and average T1 relaxation time in 14 subcortical structures
  * - Fluid-attenuated inversion recovery (FLAIR)
    - Normalized vertexwise neocortical & hippocampal FLAIR intensity, and average FLAIR intensity in 14 subcortical structures
  * - Diffusion-weighted imaging (DWI)
    - Vertexwise neocortical & hippocampal mean diffusivity and fractional anisotropy, and average mean diffusivity and fractional anisotropy in 14 subcortical structures


Data format
--------------------------------------------------------

``z-brains`` requires that your data be formatted in accordance with the `Brain Imaging Data Structure (BIDS) <https://bids.neuroimaging.io>`_ format. 
You can find the information about the BIDS specification `here <https://bids-specification.readthedocs.io/en/stable/>`_. 
We strongly recommend that you validate your data structure after the conversion, notably using the `BIDS Validator <https://bids-standard.github.io/bids-validator/>`_.


Usage overview
--------------------------------------------------------

To see a list of available options type: ::

    z-brains -help

The following help menu should then pop up with a description of available arguments for a customized run of ``z-brains`` :

  .. figure:: ../../figures/help_menu.png
	:height: 760
	:width: 700

Additional optional arguments are also available: 

  .. figure:: ../../figures/help_menu2.png
	:height: 115
	:width: 700

Code snippets relevant to running each module can be found in the appropriate section in the menu. 

See `Pre-processing requirements` for requirements on data formatting and minimal pre-procecessing to run ``z-brains``. 
Note that pre-processing relies heavily on `micapipe <https://https://micapipe.readthedocs.io/>`_ and `hippunfold <https://https://hippunfold.readthedocs.io/>`_.

The section `Feature post-processing` includes code snippets to map and smooth whole-brain feature maps along neocortical and hippocampal surface meshes, as well as sample mean volume and intensity values within each subcortical structure. 

We provide code examples to perform feature normalization in the `Feature normalization` section, including z-scoring and normative modelling approaches.

Lastly, instructions on how to generate summary reports are presented in the `Summary reports` section. 
