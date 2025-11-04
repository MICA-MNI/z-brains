"""
Z-Brains: A neuroimaging processing and analysis pipeline.
Developed by the MICA lab at McGill University.
"""

from .dataset import zbdataset, demographics
from .environment import zbenv
from .constants import Analysis, Approach, Resolution, Structure, Feature
from .clinical_reports import generate_clinical_report

__all__ = [
    "zbdataset",
    "demographics",
    "zbenv",
    "Analysis",
    "Approach",
    "Resolution",
    "Structure",
    "Feature",
    "generate_clinical_report",
]

__version__ = "2.0.0"