# 3_annotation_networks/1_map_roi_to_HO.py

import json
from pathlib import Path

import numpy as np
import pandas as pd

from config import LABELS_CSV, stage_artifacts

# ------------------------ Artifacts & Folders ------------------------
ART_CONN = stage_artifacts("2_connectivity_extraction")
CONN_DIR = ART_CONN / "connectivity"
ART_ANN = stage_artifacts("3_annotation_networks")
ART_ANN.mkdir(parents=True, exist_ok=True)


# ------------------------ Helpers ------------------------
def extract_best_label(raw_str: str) -> str:
    """
    Parse a Harvard–Oxford label cell and return the top-scoring label.

    The expected raw string format may look like:
        ["Frontal Pole": 0.65]["Insular Cortex": 0.20] ...

    If parsing fails, returns a trimmed fallback; if empty/NA, returns "Unknown".

    Parameters
    ----------
    raw_str : str
        Raw label string.

    Returns
    -------
    str
        Best label (highest score) or a fallback.
    """
    import re

    if pd.isna(raw_str) or not isinstance(raw_str, str):
        return "Unknown"

    matches = re.findall(r'\["(.+?)":\s*([\d.]+)\]', raw_str)
    if not matches:
        s = raw_str.strip()
        return s if s else "Unknown"

    lab, _ = max(((lab.strip(), float(sc)) for lab, sc in matches), key=lambda x: x[1])
    return lab


# ------------------------ Main ------------------------
def main():
    """
    Build a ROI → Harvard–Oxford label mapping aligned to the connectivity cube.

    Reads the connectivity shape to get n_rois, then aligns the provided
    ROI labels CSV (expects columns 'ROI_number' and 'Harvard-Oxford'),
    extracting the highest-probability HO label for each ROI.

    Outputs
    -------
    artifacts/3_annotation_networks/roi_to_HO.csv
        CSV with columns: ROI_number (1..n_rois), HO_label
    """
    # ------------------------ Read n_rois from connectivity metadata ------------------------
    shape = json.loads((CONN_DIR / "conn_shape.json").read_text())["shape"]
    n_rois = int(shape[1])

    # ------------------------ Load labels CSV ------------------------
    labels_path = Path(LABELS_CSV)
    assert labels_path.exists(), f"Not found: {labels_path}"
    df = pd.read_csv(labels_path)
    df.columns = df.columns.str.strip()
    assert "ROI_number" in df.columns and "Harvard-Oxford" in df.columns, \
        "CSV must contain 'ROI_number' and 'Harvard-Oxford' columns."

    df["ROI_number"] = df["ROI_number"].astype(int)
    df["HO_label"] = df["Harvard-Oxford"].apply(extract_best_label)

    # ------------------------ Align to 1..n_rois and fill missing ------------------------
    out = pd.DataFrame({"ROI_number": np.arange(1, n_rois + 1, dtype=int)})
    out = out.merge(df[["ROI_number", "HO_label"]], on="ROI_number", how="left")
    out["HO_label"] = out["HO_label"].fillna("Unknown")

    # ------------------------ Save mapping ------------------------
    out_path = ART_ANN / "roi_to_HO.csv"
    out.to_csv(out_path, index=False)
    print(f"✓ ROI → Harvard–Oxford mapping written → {out_path}")


# ------------------------ Entrypoint ------------------------
if __name__ == "__main__":
    main()
