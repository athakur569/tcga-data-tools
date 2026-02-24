"""Group TCGA slides by patient (case) for patient-level splitting.

TCGA patients often have multiple slides (e.g. primary tumor + metastasis, or
repeat sampling). Slide-level splitting on multi-slide patients causes data
leakage: the same patient appears in both train and test sets.

PatientGrouper collapses slide-level metadata to patient-level, enforcing
the invariant that each patient appears exactly once.
"""

from __future__ import annotations

import logging
import warnings

import pandas as pd

logger = logging.getLogger(__name__)

_PATIENT_PATTERN = r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})"


class RepeatedSlidesWarning(UserWarning):
    """Raised when a patient has more than one slide in the same split cohort."""


class PatientGrouper:
    """Collapse slide-level TCGA data to patient-level.

    Args:
        patient_col:  Column name for the patient/case identifier.
        slide_col:    Column name for the slide identifier.
        label_col:    Column name for the class label (for multi-label checks).
    """

    def __init__(
        self,
        patient_col: str = "case_id",
        slide_col: str = "slide_id",
        label_col: str | None = "label",
    ) -> None:
        self.patient_col = patient_col
        self.slide_col = slide_col
        self.label_col = label_col

    def infer_patient_id(self, df: pd.DataFrame, id_col: str = "submitter_id") -> pd.DataFrame:
        """Infer patient ID from TCGA submitter IDs (first 12 characters).

        TCGA-XX-XXXX is the patient identifier; everything after is sample/analyte info.
        """
        df = df.copy()
        df[self.patient_col] = df[id_col].str.extract(_PATIENT_PATTERN, expand=False)
        missing = df[self.patient_col].isna().sum()
        if missing > 0:
            logger.warning("%d rows could not have patient ID inferred", missing)
        return df

    def collapse(
        self,
        df: pd.DataFrame,
        agg: dict | None = None,
    ) -> pd.DataFrame:
        """Collapse to one row per patient.

        Args:
            df:  Slide-level DataFrame with patient_col already populated.
            agg: Custom aggregation dict passed to groupby.agg().
                 Default: list all slide IDs, take first label, sum file sizes.

        Returns:
            Patient-level DataFrame.
        """
        if self.patient_col not in df.columns:
            raise ValueError(
                f"Column '{self.patient_col}' not found. "
                "Call infer_patient_id() first or set patient_col correctly."
            )

        if agg is None:
            agg = {self.slide_col: list}
            if self.label_col and self.label_col in df.columns:
                agg[self.label_col] = "first"
            if "file_size" in df.columns:
                agg["file_size"] = "sum"
            if "tss" in df.columns:
                agg["tss"] = "first"

        patient_df = df.groupby(self.patient_col).agg(agg).reset_index()
        patient_df["n_slides"] = patient_df[self.slide_col].apply(len)

        multi_slide = (patient_df["n_slides"] > 1).sum()
        if multi_slide > 0:
            warnings.warn(
                f"{multi_slide} patients have multiple slides. "
                "Ensure splits are patient-grouped, not slide-grouped.",
                RepeatedSlidesWarning,
                stacklevel=2,
            )

        logger.info(
            "Collapsed %d slides → %d patients (%d multi-slide)",
            len(df), len(patient_df), multi_slide,
        )
        return patient_df

    def verify_no_leakage(
        self,
        train_ids: "list | set",
        val_ids: "list | set",
        test_ids: "list | set",
    ) -> bool:
        """Assert no patient appears in more than one split."""
        train, val, test = set(train_ids), set(val_ids), set(test_ids)
        tv = train & val
        tt = train & test
        vt = val & test
        if tv or tt or vt:
            raise ValueError(
                f"Split leakage detected! "
                f"train∩val={len(tv)}, train∩test={len(tt)}, val∩test={len(vt)}"
            )
        return True
