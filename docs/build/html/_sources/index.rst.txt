.. MICAPIPE documentation master file, created by
   sphinx-quickstart on Wed Jul 15 16:09:38 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

   <style type="text/css">
      hr {
      width: 100%;
      height: 1px;
      background-color: #e8e8e8;
      margin-top: 24px;
      }
   </style>

.. .. image:: https://readthedocs.org/projects/z-brains/badge/?version=latest
  :target: https://z-brains.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status

**Welcome to z-brains's documentation!**
========================================================

`z-brains` identifies patient-specific signal alterations in brain morphology and microstructure, 
using features with previously demonstrated potential for localization and lateralization of epileptogenic lesions.

.. raw:: html

   <br>

___________________________________________________________________________________________________


Breaking news ️‍📰
--------------------------------------------------------
`z-brains version 0.0.2 "Reborn"` is up and running!


About 👁️‍🗨️
--------------------------------------------------------
`z-brains` is an open toolbox promoting standardized procedures to map and normalize structural imaging features in epilepsy research.
This pipeline integrates processing streams to map and normalize T1-weighted (cortical thickness, subcortical volume), diffusion-weighted (mean diffusivity, fractional anisotropy),
FLAIR (intensity), and quantitative T1 (intensity) imaging data to standard neocortical, subcortical, and hippocampal spaces.

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting started

   pages/01.install/index
   pages/02.tutorialinstructions/index
   pages/03.whatsnew/index
   
.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Pre-processing requirements

   pages/04.cortexstcx/index
   pages/05.hippocampus/index
   
.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Feature post-processing

   pages/06.cortex/index
   pages/07.subcortex/index
   pages/08.hippocampus/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Feature normalization

   pages/09.zscoring/index
   pages/10.normative/index
   
.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Summary reports

   pages/11.features/index
   pages/12.reports/index
   

___________________________________________________________________________________________________

.. raw:: html

   <br>

Core development team 🧠
--------------------------------------------------------

`z-brains` is developed by members of the MICA-lab (https://mica-mni.github.io) and collaborators 
at the McConnell Brain Imaging Centre of the Montreal Neurological Institute and Harvard Medical School.

- **Jessica Royer**, *MICA Lab - Montreal Neurological Institute*
- **Sara Larivière**, *Brigham and Women's Hospital | Harvard Medical School*
- **Raúl Rodríguez-Cruces**, *MICA Lab - Montreal Neurological Institute*
- **Oualid Benkarim**, *MICA Lab - Montreal Neurological Institute*
- **Jordan Dekraker**, *MICA Lab - Montreal Neurological Institute*
- **Judy Chen**, *MICA Lab - Montreal Neurological Institute*
- **Boris Bernhardt**, *MICA Lab - Montreal Neurological Institute*
