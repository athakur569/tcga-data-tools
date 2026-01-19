"""GDC API client for fetching TCGA WSI manifests and case metadata."""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)

_GDC_API = "https://api.gdc.cancer.gov"
_FILES_ENDPOINT = f"{_GDC_API}/files"
_CASES_ENDPOINT = f"{_GDC_API}/cases"


class GDCClient:
    """Query the GDC API for TCGA whole-slide image manifests.

    Args:
        timeout:     HTTP request timeout in seconds.
        retry_delay: Seconds to wait between retries on failure.
        max_retries: Max number of retries per request.
    """

    def __init__(
        self,
        timeout: int = 30,
        retry_delay: float = 2.0,
        max_retries: int = 3,
    ) -> None:
        self.timeout = timeout
        self.retry_delay = retry_delay
        self.max_retries = max_retries

    def get_manifest(
        self,
        project: str,
        data_type: str = "Slide Image",
        experimental_strategy: str = "Tissue Slide",
        page_size: int = 500,
    ) -> pd.DataFrame:
        """Fetch WSI file manifest for a TCGA project.

        Args:
            project:                 e.g. 'TCGA-BRCA', 'TCGA-LUAD', 'TCGA-COAD'
            data_type:               GDC data type filter.
            experimental_strategy:   GDC experimental strategy filter.
            page_size:               Records per API page.

        Returns:
            DataFrame with columns: file_id, file_name, case_id, slide_id,
            submitter_id, md5sum, file_size.
        """
        filters = {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
                {"op": "=", "content": {"field": "data_type", "value": data_type}},
                {"op": "=", "content": {"field": "experimental_strategy", "value": experimental_strategy}},
            ],
        }
        fields = [
            "file_id", "file_name", "md5sum", "file_size",
            "cases.case_id", "cases.submitter_id",
            "cases.samples.portions.slides.submitter_id",
        ]

        rows, from_ = [], 0
        total = None

        with tqdm(desc=f"GDC manifest {project}", unit="files") as pbar:
            while True:
                data = self._request(
                    _FILES_ENDPOINT,
                    params={
                        "filters": __import__("json").dumps(filters),
                        "fields": ",".join(fields),
                        "format": "JSON",
                        "size": page_size,
                        "from": from_,
                    },
                )
                hits = data.get("data", {}).get("hits", [])
                if total is None:
                    total = data.get("data", {}).get("pagination", {}).get("total", 0)
                    pbar.total = total

                for hit in hits:
                    case = hit.get("cases", [{}])[0]
                    rows.append({
                        "file_id": hit.get("file_id", ""),
                        "file_name": hit.get("file_name", ""),
                        "md5sum": hit.get("md5sum", ""),
                        "file_size": hit.get("file_size", 0),
                        "case_id": case.get("case_id", ""),
                        "submitter_id": case.get("submitter_id", ""),
                    })
                pbar.update(len(hits))
                from_ += page_size
                if from_ >= (total or 0):
                    break

        logger.info("Fetched %d files for %s", len(rows), project)
        return pd.DataFrame(rows)

    def _request(self, url: str, params: dict, attempt: int = 0) -> dict[str, Any]:
        try:
            r = requests.get(url, params=params, timeout=self.timeout)
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, ValueError) as e:
            if attempt < self.max_retries:
                logger.warning("Request failed (%s), retrying in %.1fs...", e, self.retry_delay)
                time.sleep(self.retry_delay)
                return self._request(url, params, attempt + 1)
            raise
