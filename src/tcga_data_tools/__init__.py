"""tcga-data-tools: TCGA data access and preprocessing utilities."""

from .gdc_client import GDCClient
from .tss_parser import TSSParser
from .patient_grouper import PatientGrouper
from .cohort_stats import cohort_summary

__version__ = "0.1.0"
__all__ = ["GDCClient", "TSSParser", "PatientGrouper", "cohort_summary"]
