# 2_connectivity_extraction/1_compute_connectivity.py

import numpy as np
import pandas as pd
from pathlib import Path
from config import DATA_DIR, stage_artifacts, MAX_BAD_ROIS_RATIO
from utils_io import read_1d, save_memmap_cube, dump_json
from utils_conn import corr_from_timeseries

# ------------------------ Output folders ------------------------
ART = stage_artifacts("2_connectivity_extraction")
PHENO_ART = stage_artifacts("1_data_cleaning")   # cleaned phenotype
OUT_DIR = ART / "connectivity"                   # main outputs
QC_DIR = ART / "quality_check"                   # per-subject QC
OUT_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------ Helpers ------------------------
def fisher_z(R, eps: float = 1e-6) -> np.ndarray:
    """
    Apply Fisher z-transform to a correlation matrix.

    Parameters
    ----------
    R : np.ndarray
        Correlation matrix with values in [-1, 1].
    eps : float, optional
        Small epsilon to clip away from ±1 to avoid infinities.

    Returns
    -------
    np.ndarray
        Fisher z-transformed matrix.
    """
    R = np.clip(R, -1 + eps, 1 - eps)  # avoid ±inf
    return 0.5 * np.log((1.0 + R) / (1.0 - R))


# ------------------------ Main ------------------------
def main():
    """
    Build subject-level connectivity matrices (Fisher z) from ROI time series.

    Workflow
    --------
    1) Read cleaned phenotype and subject list.
    2) Infer n_rois from the first readable .1D file.
    3) For each subject:
       - Read ROI time series (.1D).
       - Compute correlation matrix with `corr_from_timeseries`.
       - Reject subjects with too many constant ROIs (> MAX_BAD_ROIS_RATIO).
       - Store Fisher z matrix, label (DX_GROUP) and QC info.
    4) Save results:
       - Memmap cube `conn.memmap` of shape (N, n_rois, n_rois).
       - `subjects.npy`, `labels.npy`, `roi_ids.npy`, `n_rois.txt`.
       - QC CSV and a small JSON with shape metadata.
    """
    # ------------------------ Load cleaned phenotype ------------------------
    subjects = pd.read_csv(PHENO_ART / "subjects.csv")["FILE_ID"].astype(str).tolist()
    pheno = pd.read_parquet(PHENO_ART / "pheno.parquet")
    assert pheno.index.name == "FILE_ID", "pheno.parquet must be indexed by FILE_ID."

    # ------------------------ Determine n_rois ------------------------
    n_rois = None
    for sid in subjects:
        f = Path(DATA_DIR) / f"{sid}_rois_cc400.1D"
        if f.exists():
            try:
                n_rois = read_1d(f).shape[1]
                break
            except Exception:
                continue
    assert n_rois is not None, "Unable to determine n_rois."
    print(f"n_rois = {n_rois}")

    # ------------------------ Per-subject loop ------------------------
    corr_list, ok_subjects, labels, qc_rows = [], [], [], []
    for sid in subjects:
        f = Path(DATA_DIR) / f"{sid}_rois_cc400.1D"
        if not f.exists():
            qc_rows.append({"FILE_ID": sid, "status": "missing_file"})
            continue

        try:
            ts = read_1d(f).values
        except Exception:
            qc_rows.append({"FILE_ID": sid, "status": "read_error"})
            continue

        if ts.ndim != 2 or ts.shape[1] != n_rois:
            qc_rows.append({"FILE_ID": sid, "status": "shape_mismatch"})
            continue

        corr, rep = corr_from_timeseries(ts, report=True)
        n_bad = int(rep.get("n_bad_rois", 0))
        bad_ratio = n_bad / float(n_rois)
        if bad_ratio > MAX_BAD_ROIS_RATIO:
            print(f"[skip] too many constant ROIs {sid}: {n_bad}/{n_rois} ({bad_ratio:.1%})")
            qc_rows.append({"FILE_ID": sid, "status": "too_many_bad_rois", "n_bad": n_bad})
            continue

        M = fisher_z(corr)  # Fisher z-transform
        corr_list.append(M.astype("float32"))
        ok_subjects.append(sid)
        labels.append(int(pheno.loc[sid, "DX_GROUP"]))
        qc_rows.append({"FILE_ID": sid, "status": "kept", "n_bad": n_bad})

    print(f"✓ kept subjects: {len(ok_subjects)}/{len(subjects)}")

    # ------------------------ Save memmap cube ------------------------
    conn_mm_path = OUT_DIR / "conn.memmap"
    conn = save_memmap_cube(conn_mm_path, shape=(len(ok_subjects), n_rois, n_rois))
    for i, C in enumerate(corr_list):
        conn[i] = C

    # ------------------------ QC & metadata ------------------------
    pd.DataFrame(qc_rows).to_csv(QC_DIR / "qc_connectivity_per_subject.csv", index=False)
    print("↳ detailed QC written -> artifacts/2_connectivity_extraction/quality_check/qc_connectivity_per_subject.csv")
    np.save(OUT_DIR / "subjects.npy", np.array(ok_subjects, dtype=object))
    np.save(OUT_DIR / "labels.npy", np.array(labels, dtype=int))
    np.save(OUT_DIR / "roi_ids.npy", np.arange(n_rois, dtype=int))
    (OUT_DIR / "n_rois.txt").write_text(str(n_rois))
    dump_json(OUT_DIR / "conn_shape.json", {"shape": [len(ok_subjects), n_rois, n_rois]})

    print(f"✓ conn.memmap written ({len(ok_subjects)} subjects, {n_rois} ROIs)")


# ------------------------ Entrypoint ------------------------
if __name__ == "__main__":
    main()
