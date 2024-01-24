.. figure:: ./data/zbrains_banner.png
   :alt: alternate text
   :align: center

Multimodal lesion mapping in focal epilepsy with ``zbrains``
--------------------------------------------

.. image:: https://img.shields.io/github/v/tag/MICA-MNI/z-brains
  :target: https://github.com/MICA-MNI/z-brains
  :alt: version

.. image:: https://img.shields.io/badge/license-BSD-brightgreen
   :target: https://opensource.org/licenses/BSD-3-Clause

.. image:: https://readthedocs.org/projects/z-brains/badge/?version=latest&color=brightgreen
  :target: https://z-brains.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status

.. image:: https://img.shields.io/github/issues/MICA-MNI/z-brains?color=brightgreen
  :target: https://github.com/MICA-MNI/z-brains/issues
  :alt: GitHub issues

.. image:: https://img.shields.io/github/stars/MICA-MNI/z-brains.svg?style=flat&label=%E2%9C%A8%EF%B8%8F%20be%20a%20stargazer&color=brightgreen
    :target: https://github.com/MICA-MNI/z-brains/stargazers
    :alt: GitHub stars

\

``zbrains`` is developed by `MICA-Lab <https://mica-mni.github.io>`_, affiliated with `McGill University <https://www.mcgill.ca/>`_, the Montreal Neurological Institute and Hospital "`the Neuro <https://www.mcgill.ca/neuro/>`_," and the McConnell Brain Imaging Center (`BIC <https://www.mcgill.ca/bic/>`_).

   This open access processing and analysis tool aims identify patient-specific anomalies in brain morphology and microstructure, using features with previously demonstrated potential to accurately localize epileptogenic lesions.
   ``zbrains`` uses a set of known software dependencies developed by other groups and aggregated in a published pipeline `micapipe <https://github.com/MICA-MNI/micapipe>`_.

.. Installation
.. --------------------------------------------

.. Make sure set MICAPIPE and ZBRAINS variables, and add their function to your PATH. For example:
.. .. code-block bash::
..    export MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
..    export PATH=${PATH}:${MICAPIPE}:${MICAPIPE}/functions
..    source ${MICAPIPE}/functions/init.sh

..    export ZBRAINS=/data/mica1/03_projects/jordand/z-brains
..    export PATH=${PATH}:${ZBRAINS}:${ZBRAINS}/functions

.. ::

Tutorial
--------------------------------------------

You must be inside the ``zbrains`` repository to run the following commands.

.. code-block:: bash

   cd /path/to/zbrains/repository

``zbrains`` requires input and output directories:

- ``pth_dataset`` points to the BIDS-format dataset path where micapipe, hippunfold, and zbrains directories will be stored
- ``micapipe_dir`` contains the output of ``micapipe`` previously run on the BIDS dataset (also BIDS-format)
- ``hippunfold_dir`` contains the output of ``hippunfold`` previously run on the BIDS dataset (also BIDS-format)
- ``zbrains_dir`` points to the directory that will hold ``z-brains`` outputs

.. code-block:: bash

    # path for dataset in BIDS structure
    pth_dataset=/path/to/BIDS_dataset

    micapipe_dir="name_of_micapipe_folder"
    hippunfold_dir="name_of_hippunfold_folder"
    zbrains_dir="name_of_z-brains_folder"


Preparing control data
---------------------------------------------

A ``.csv`` file with ID and session for healthy controls are required

.. code-block:: bash

  # csv file with ID and session for control participants to be processed
  demo_controls='/path/to/control/participants.csv'

To process features for healthy controls as a batch, run the following. Note that column names of the healthy controls ``.csv`` file must be specified under ``--column_map``.

.. code-block:: bash
  
      ./zbrains_batch --run proc \
         --demo "${demo_controls}" \
         --dataset "${pth_dataset}" \
         --zbrains ${zbrains_dir} \
         --micapipe ${micapipe_dir} \
         --hippunfold ${hippunfold_dir} \
         --column_map participant-id=ID session_id=SES \
         --verbose 2 \
         --scheduler_options "-n 20" #specify threads here


Processing and analyzing patient features
------------------------------------------------

.. code-block:: bash

    # specify the list of subject IDs along with corresponding session
    subject_ids=(sub-PX001 sub-PX002)
    session_ids=(ses-01 ses-01)

    # csv file with ID and session for control participants for comparison
    demo_controls='/path/to/control/participants.csv'

      for i in "${!subject_ids[@]}"
      do
        sid=${subject_ids[$i]}
        ses=${session_ids[$i]:-""}
        echo -e "$i\t\tID: $sid, SES: $ses"
      
      
        ./zbrains --run analysis \
                  --sub "${sid}" \
                  --ses "${ses}" \
                  --dataset ${pth_dataset} \
                  --zbrains ${zbrains_dir} \
                  --demo_ref ${demo_controls} \
                  --dataset_ref ${pth_dataset} \
                  --zbrains_ref ${zbrains_dir} \
                  --column_map participant_id=ID session_id=SES \
                  --verbose 2
      done
