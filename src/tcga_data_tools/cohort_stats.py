"""Cohort summary statistics for TCGA manifests."""

from __future__ import annotations

import pandas as pd


def cohort_summary(
    df: pd.DataFrame,
    patient_col: str = "case_id",
    label_col: str | None = "label",
    tss_col: str | None = "tss",
    file_size_col: str | None = "file_size",
) -> dict:
    """Compute summary statistics for a TCGA cohort DataFrame.

    Returns a dict with:
    - n_slides, n_patients
    - label distribution (if label_col present)
    - institution distribution (if tss_col present)
    - total_gb (if file_size_col present)
    - multi_slide_patients count
    """
    stats: dict = {
        "n_slides": len(df),
        "n_patients": df[patient_col].nunique() if patient_col in df.columns else None,
    }

    if label_col and label_col in df.columns:
        stats["label_counts"] = df[label_col].value_counts().to_dict()
        stats["label_balance"] = (
            df.groupby(patient_col)[label_col].first().value_counts().to_dict()
            if patient_col in df.columns else {}
        )

    if tss_col and tss_col in df.columns:
        stats["n_institutions"] = df[tss_col].nunique()
        stats["institution_counts"] = df[tss_col].value_counts().to_dict()

    if file_size_col and file_size_col in df.columns:
        stats["total_gb"] = round(df[file_size_col].sum() / 1e9, 2)

    if patient_col in df.columns and "slide_id" in df.columns:
        slides_per_patient = df.groupby(patient_col)["slide_id"].count()
        stats["multi_slide_patients"] = int((slides_per_patient > 1).sum())
        stats["max_slides_per_patient"] = int(slides_per_patient.max())

    return stats


def print_cohort_summary(df: pd.DataFrame, cohort_name: str = "Cohort", **kwargs) -> None:
    stats = cohort_summary(df, **kwargs)
    print(f"\n{'='*50}")
    print(f"  {cohort_name}")
    print(f"{'='*50}")
    print(f"  Slides    : {stats['n_slides']}")
    if stats.get("n_patients"):
        print(f"  Patients  : {stats['n_patients']}")
    if "label_counts" in stats:
        print(f"  Labels    : {stats['label_counts']}")
    if "n_institutions" in stats:
        print(f"  Sites     : {stats['n_institutions']} TSS codes")
    if "total_gb" in stats:
        print(f"  Total size: {stats['total_gb']} GB")
    if "multi_slide_patients" in stats:
        print(f"  Multi-slide patients: {stats['multi_slide_patients']}")
    print(f"{'='*50}\n")
