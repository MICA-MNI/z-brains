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

    
Installation
--------------------------------------------

Make sure set MICAPIPE and ZBRAINS variables, and add their function to your PATH. For example:
```
export MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
export PATH=${PATH}:${MICAPIPE}:${MICAPIPE}/functions
source ${MICAPIPE}/functions/init.sh

export ZBRAINS=/data/mica1/03_projects/jordand/z-brains
export PATH=${PATH}:${ZBRAINS}:${ZBRAINS}/functions
```
