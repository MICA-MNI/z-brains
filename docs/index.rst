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

**Welcome to zbrains's documentation!**
========================================================

\

.. image:: https://img.shields.io/badge/license-BSD-brightgreen
   :target: https://opensource.org/licenses/BSD-3-Clause

.. image:: https://img.shields.io/github/v/tag/MICA-MNI/z-brains
  :target: https://github.com/MICA-MNI/z-brains
  :alt: version

.. image:: https://readthedocs.org/projects/z-brains/badge/?version=latest&color=brightgreen
  :target: https://z-brains.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status

.. image:: https://img.shields.io/github/issues/MICA-MNI/z-brains?color=brightgreen
  :target: https://github.com/MICA-MNI/z-brains/issues
  :alt: GitHub issues

.. image:: https://img.shields.io/github/stars/MICA-MNI/z-brains.svg?style=flat&label=%E2%9C%A8%EF%B8%8F%20be%20a%20stargazer&color=brightgreen
    :target: https://github.com/MICA-MNI/z-brains/stargazers
    :alt: GitHub stars

`zbrains` identifies patient-specific signal alterations in brain morphology and microstructure,
using features with previously demonstrated potential for localization and lateralization of epileptogenic lesions.

.. raw:: html

   <br>


Breaking news ️‍📰
--------------------------------------------------------
``zbrains`` version v1.0 "lido" is up and running!


About 👁️‍🗨️
--------------------------------------------------------
``zbrains`` is an open toolbox promoting standardized procedures to map and normalize structural imaging features in epilepsy research.
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

   pages/04.micapipe/index
   pages/05.hippunfold/index

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
   :caption: Main outputs

   pages/11.features/index
   pages/12.reports/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Adding plugin features

   pages/13.plugins/index



___________________________________________________________________________________________________

.. raw:: html

   <br>

Core development team 🧠
--------------------------------------------------------

`z-brains` is developed by members of the MICA-lab (https://mica-mni.github.io) and collaborators
at the McConnell Brain Imaging Centre of the Montreal Neurological Institute and Harvard Medical School.

- **Ian Goodall-Halliwell**, *MICA Lab - Montreal Neurological Institute*
- **Jessica Royer**, *MICA Lab - Montreal Neurological Institute*
- **Sara Larivière**, *Brigham and Women's Hospital | Harvard Medical School*
- **Raúl Rodríguez-Cruces**, *MICA Lab - Montreal Neurological Institute*
- **Oualid Benkarim**, *MICA Lab - Montreal Neurological Institute*
- **Jordan Dekraker**, *MICA Lab - Montreal Neurological Institute*
- **Yigu Zhou**, *MICA Lab - Montreal Neurological Institute*
- **Youngeun Hwang**, *MICA Lab - Montreal Neurological Institute*
- **Judy Chen**, *MICA Lab - Montreal Neurological Institute*
- **Boris Bernhardt**, *MICA Lab - Montreal Neurological Institute*
