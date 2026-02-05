"""Extract TSS (Tissue Source Site) codes from TCGA slide IDs.

TCGA slide IDs follow the format:
    TCGA-{TSS}-{participant}-{sample}{vial}-{portion}{analyte}-{plate}-{center}

The TSS is the 2-character code identifying the submitting institution.
Extracting it enables institution-stratified splitting to prevent site leakage.

Example:
    TCGA-E2-A1B5-01Z-00-DX1  →  TSS = 'E2'
    TCGA-A8-A075-01A-01-BS1  →  TSS = 'A8'
"""

from __future__ import annotations

import re

import pandas as pd

_TSS_PATTERN = re.compile(r"^TCGA-([A-Z0-9]{2})-", re.IGNORECASE)

# Partial lookup: TSS → institution name (subset; extend as needed)
TSS_LOOKUP: dict[str, str] = {
    "E2": "International Genomics Consortium",
    "A1": "University of Alabama at Birmingham",
    "A2": "Asterand Bioscience",
    "A7": "Cureline",
    "A8": "Indivumed",
    "AC": "International Genomics Consortium",
    "AN": "Asterand Bioscience",
    "AO": "Cureline",
    "AQ": "BioreclamationIVT",
    "AR": "Cureline",
    "B6": "Roswell Park Cancer Institute",
    "BH": "BioreclamationIVT",
    "D8": "ILSbio",
    "E9": "Nationwide Children's Hospital",
    "EW": "Cureline",
    "GI": "Tissue Solutions",
    "GM": "OSU",
    "OL": "Asterand Bioscience",
    "PL": "Mayo Clinic",
}


class TSSParser:
    """Extract TSS codes from TCGA slide/case IDs."""

    def extract_one(self, slide_id: str) -> str | None:
        """Extract TSS from a single TCGA slide ID.

        Returns:
            2-character TSS string (uppercase), or None if not a valid TCGA ID.
        """
        m = _TSS_PATTERN.match(str(slide_id))
        return m.group(1).upper() if m else None

    def extract(self, slide_ids: "pd.Series | list[str]") -> pd.Series:
        """Extract TSS from a Series or list of slide IDs.

        Returns:
            pd.Series of TSS strings (NaN where extraction fails).
        """
        if not isinstance(slide_ids, pd.Series):
            slide_ids = pd.Series(slide_ids)
        return slide_ids.str.extract(r"^TCGA-([A-Z0-9]{2})-", expand=False).str.upper()

    def add_institution(
        self,
        df: pd.DataFrame,
        tss_col: str = "tss",
        out_col: str = "institution",
    ) -> pd.DataFrame:
        """Map TSS codes to institution names using the built-in lookup."""
        df = df.copy()
        df[out_col] = df[tss_col].map(TSS_LOOKUP).fillna("Unknown")
        return df

    def validate(self, df: pd.DataFrame, id_col: str = "slide_id") -> dict:
        """Report how many rows have valid vs missing TSS codes."""
        tss = self.extract(df[id_col])
        valid = tss.notna().sum()
        return {
            "total": len(df),
            "valid_tss": int(valid),
            "missing_tss": int(len(df) - valid),
            "unique_tss": int(tss.nunique()),
            "tss_counts": tss.value_counts().to_dict(),
        }
