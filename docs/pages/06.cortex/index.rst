.. _postcortex:

.. title:: Cortical post-processing

Cortical feature post-processing
============================================================

This section describes how to generate raw and normalized features sampled along the neocortical surface generated by FreeSurfer or FastSurfer via micapipe.


Define your directories: Cortex
--------------------------------------------------------

zbrains requires input and output directories:

``root_path`` points to the BIDS-format dataset that stores imaging data
``rawdir`` contains the raw imaging data
``micapipedir`` contains the output of micapipe previously run on the BIDS dataset
``hippdir`` contains the output of hippunfold previously run on the BIDS dataset
``outdir`` points to the directory that will hold zbrains outputs

In practice, these should be declared as variables:

    .. code-block:: bash
        :caption: Declaring path variables
        :linenos:

        # path for dataset in BIDS structure
        root_path=/path/to/BIDS_dataset
        rawdir=${root_path}/rawdata
        micapipedir=${root_path}/derivatives/micapipe_folder
        hippdir=${root_path}/derivatives/hippunfold_folder
        outdir=${root_path}/derivatives/zbrains_folder


Run cortical post-processing
--------------------------------------------------------

After defining our directories, we can move on to processing our first patient.
Let's declare which patient we want to process (with their ID in the BIDS directory, omitting the 'sub-' and the session that should be processed if necessary).

    .. code-block:: bash
        :caption: Declaring subject variables
        :linenos:

        id=PX001
        ses=01

To process the subject, we can specify the command as follows:

    .. code-block:: bash
        :caption: Basic zbrains run: cortical processing
        :linenos:

        zbrains -sub "$id" -ses "$ses" \
            -rawdir "${rawdir}" \
            -micapipedir "${micapipedir}" \
            -hippdir "${hippdir}" \
            -outdir "${outdir}" \
            -run proc \
            -struct cortex \
            -verbose 2

This will generate surface-based cortical features at the default smoothing (5mm) and mesh resolution (both 5k and 32k) for all available modalities.

.. admonition:: Customize your zbrains run!

	A list of options and flags can be specified for a more personalized run of zbrains. Check out the help menu ``zbrains -help``
