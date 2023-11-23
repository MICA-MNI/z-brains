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

    
Multimodal lesion mapping in focal epilepsy with `z-brains`
--------------------------------------------

`z-brains` is developed by `MICA-Lab <https://mica-mni.github.io>`_, affiliated with `McGill University <https://www.mcgill.ca/>`_, the Montreal Neurological Institute and Hospital "`the Neuro <https://www.mcgill.ca/neuro/>`_," and the McConnell Brain Imaging Center (`BIC <https://www.mcgill.ca/bic/>`_).

This open access processing and analysis tool aims identify patient-specific anomalies in brain morphology and microstructure, using features with previously demonstrated potential to accurately localize epileptogenic lesions. 

`z-brains` uses a set of known software dependencies developped by other groups and aggregated in a published pipeline `micapipe <https://github.com/MICA-MNI/micapipe>`_.

    
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

ZBRAINS requires input and output directories:
   `root_path` points to the BIDS-format dataset that stores imaging data
   `rawdir` contains the raw imaging data
   `micapipedir` contains the output of MICAPIPE previously run on the BIDS dataset
   `hippdir` contains the output of HIPPUNFOLD previously run on the BIDS dataset
   `outdir` points to the directory that will hold ZBRAINS outputs

.. code-block:: bash

    # Path for dataset in BIDS structure
    root_path=/path/to/BIDS_dataset

    rawdir=${root_path}/rawdata
    micapipedir=${root_path}/derivatives/micapipe_folder
    hippdir=${root_path}/derivatives/hippunfold_folder
    outdir=${root_path}/derivatives/z-brains_folder


Preparing control data
---------------------------------------------

To process features for healthy controls, subject and session identifiers are required

.. code-block:: bash

  # csv file with ID and session for control participants to be processed
  PATH_CSV_CONTROLS='/path/to/control/participants.csv'

  while IFS=',' read -r id ses rest
  do
      ./z-brains -sub "$id" -ses "$ses" \
      -rawdir "${rawdir}" \
      -micapipedir "${micapipedir}" \
      -hippdir "${hippdir}" \
      -outdir "${outdir}" \
      -run proc \
      -mica \
      -verbose 2 

  done <<< "$(tail -n +2 "${PATH_CSV_CONTROLS}")"


Processing and analyzing patient features
------------------------------------------------

.. code-block:: bash

    # Specify the list of subject IDs along with corresponding session
    px_id=(PX001 PX002 PX003)
    px_ses=(1 1 1)

    # csv file with ID and session for control participants for comparison
    PATH_CSV_CONTROLS='/path/to/control/participants.csv'

    i=0
    for id in "${px_id[@]}"
    do
        ses=${px_ses[$i]}
        
        ./z-brains -sub "$id" -ses "$ses" \
        -rawdir "${rawdir}" \
        -micapipedir "${micapipedir}" \
        -hippdir "${hippdir}" \
        -outdir "${outdir}" \
        -approach "zscore" \
        -demo_cn "${PATH_CSV_CONTROLS}" \
        -mica -verbose 2

        i=$((i+1))
        
    done

