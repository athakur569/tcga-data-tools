# tcga-data-tools

Utilities for working with TCGA (The Cancer Genome Atlas) data via the GDC API: manifest management, slide download, patient-level grouping, TSS code extraction, and cohort statistics.

## Features
- Query the GDC API for WSI manifests by cohort (BRCA, COAD, LUAD, etc.)
- Extract TSS codes (institution identifiers) from TCGA slide IDs
- Group slides by patient — critical for patient-level train/val/test splitting
- Cohort summary statistics (case counts, institution breakdown, label distribution)

## Installation
```bash
pip install tcga-data-tools
```

## Quick start
```python
from tcga_data_tools import GDCClient, TSSParser, PatientGrouper

# Fetch manifest for TCGA-BRCA
client = GDCClient()
manifest = client.get_manifest(project="TCGA-BRCA", data_type="Slide Image")
manifest.to_csv("brca_manifest.csv", index=False)

# Extract TSS codes from slide IDs
parser = TSSParser()
df["tss"] = parser.extract(df["slide_id"])

# Group by patient
grouper = PatientGrouper()
patient_df = grouper.collapse(df)
```

## TSS codes
TSS (Tissue Source Site) codes are the 2-character institutional identifiers embedded in TCGA slide IDs: `TCGA-E2-A1B5-...` → TSS `E2`.
These enable institution-stratified splitting to prevent site leakage.

## License
MIT
