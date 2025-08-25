# 2_connectivity_extraction/2_regress_covariates.py

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from config import stage_artifacts

# ------------------------ Folders ------------------------
ART = stage_artifacts("2_connectivity_extraction")
PHENO_ART = stage_artifacts("1_data_cleaning")
CONN_DIR = ART / "connectivity"
REG_DIR = ART / "regressed"
REG_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------ Options ------------------------
WRITE_EDGES_LONG = False  # set True if you want to export edges_long.parquet (large file)


# ------------------------ Helpers ------------------------
def vectorize_connectivity(conn: np.ndarray):
    """
    Vectorize symmetric connectivity matrices (upper triangle only).

    Parameters
    ----------
    conn : np.ndarray
        Connectivity array of shape (n_subjects, n_rois, n_rois).

    Returns
    -------
    conn_vec : np.ndarray
        Vectorized upper triangle (n_subjects, n_edges).
    iu : tuple
        Indices (i, j) of the upper triangle.
    n_rois : int
        Number of ROIs.
    """
    n_subj, n_rois, _ = conn.shape
    iu = np.triu_indices(n_rois, k=1)
    return conn[:, iu[0], iu[1]], iu, n_rois


def make_design_matrix(pheno: pd.DataFrame, covariates: list) -> pd.DataFrame:
    """
    Build design matrix from covariates, with categorical encoding.

    Parameters
    ----------
    pheno : pd.DataFrame
        Phenotype dataframe aligned to subjects.
    covariates : list
        List of covariate column names.

    Returns
    -------
    X : pd.DataFrame
        Design matrix with intercept.
    """
    df = pheno[covariates].copy()
    for c in df.columns:
        if df[c].dtype == "O" or str(df[c].dtype) == "category":
            df[c] = df[c].astype("category")
    X = pd.get_dummies(df, drop_first=True, dummy_na=False)
    return sm.add_constant(X.astype(float))


def regress_out_covariates_fast(Y: np.ndarray, X: np.ndarray):
    """
    Regress out covariates from edge values using OLS (vectorized).

    Parameters
    ----------
    Y : np.ndarray
        Shape (n_subjects, n_edges). Edge values.
    X : np.ndarray
        Shape (n_subjects, p). Design matrix.

    Returns
    -------
    resid : np.ndarray
        Residuals, same shape as Y, NaN where rows are invalid.
    row_valid : np.ndarray
        Boolean mask of valid rows.
    """
    row_valid = np.isfinite(X).all(axis=1) & np.isfinite(Y).all(axis=1)
    Xv, Yv = X[row_valid], Y[row_valid]
    # Multivariate least squares: beta shape (p, n_edges)
    beta, *_ = np.linalg.lstsq(Xv, Yv, rcond=None)
    resid_v = Yv - Xv @ beta
    # Insert back residuals
    resid = np.full_like(Y, np.nan, dtype=float)
    resid[row_valid] = resid_v
    return resid, row_valid


# ------------------------ Main ------------------------
def main():
    """
    Regress out covariates (AGE, SEX, SITE) from subject-level connectivity
    matrices, saving residuals in both vectorized and full-matrix form.
    """
    # ------------------------ Load memmap ------------------------
    shape = tuple(json.loads((CONN_DIR / "conn_shape.json").read_text())["shape"])
    conn = np.memmap(CONN_DIR / "conn.memmap", dtype="float32", mode="r", shape=shape)
    ok_subjects = np.load(CONN_DIR / "subjects.npy", allow_pickle=True)

    # ------------------------ Vectorize ------------------------
    conn_vec, iu, n_rois = vectorize_connectivity(conn)  # (n_subj, n_edges)

    # ------------------------ Align phenotype ------------------------
    pheno = pd.read_parquet(PHENO_ART / "pheno.parquet")
    assert pheno.index.name == "FILE_ID", "pheno.parquet must be indexed by FILE_ID."
    pheno_aligned = pheno.reindex(ok_subjects)

    # ------------------------ Design matrix + regression ------------------------
    covariates = ["AGE_AT_SCAN", "SEX", "SITE_ID"]
    X = make_design_matrix(pheno_aligned, covariates).values
    residuals, row_valid = regress_out_covariates_fast(conn_vec.astype(float), X)

    subjects_valid = np.array(ok_subjects, dtype=object)[row_valid]
    pheno_valid = pheno_aligned.loc[subjects_valid]
    labels_valid = pheno_valid["DX_GROUP"].astype(int).to_numpy()  # 1=ASD, 2=HC
    group_name = np.where(labels_valid == 1, "ASD", "HC")

    # ------------------------ Save vectorized residuals ------------------------
    np.save(REG_DIR / "connectivity_residuals.npy", residuals)        # (n_subj, n_edges) with NaN outside row_valid
    np.save(REG_DIR / "subjects_valid.npy", subjects_valid)
    np.save(REG_DIR / "labels_valid.npy", labels_valid)
    np.save(REG_DIR / "triu_i.npy", iu[0])
    np.save(REG_DIR / "triu_j.npy", iu[1])
    (REG_DIR / "n_rois.txt").write_text(str(n_rois))
    print(f"✓ Residuals (vectorized) written → {REG_DIR}")

    # ------------------------ Save residual matrices ------------------------
    resid_v = residuals[row_valid]  # (n_valid, n_edges)
    n_valid = resid_v.shape[0]
    conn_resid = np.memmap(REG_DIR / "conn_resid.memmap", dtype="float32", mode="w+",
                           shape=(n_valid, n_rois, n_rois))
    for k in range(n_valid):
        M = np.zeros((n_rois, n_rois), dtype=np.float32)
        M[iu] = resid_v[k]
        M[(iu[1], iu[0])] = M[iu]
        np.fill_diagonal(M, 0.0)
        conn_resid[k] = M
    del conn_resid
    print(f"✓ Residual matrices written → {REG_DIR}/conn_resid.memmap ({n_valid} subjects)")

    # ------------------------ Optional long format ------------------------
    if WRITE_EDGES_LONG:
        rows = []
        for k, sid in enumerate(subjects_valid):
            v = resid_v[k]
            dfk = pd.DataFrame({
                "subject": sid,
                "group": group_name[k],
                "roi_i": iu[0].astype(np.int32),
                "roi_j": iu[1].astype(np.int32),
                "value": v.astype(np.float32),
            })
            rows.append(dfk)
        edges_long = pd.concat(rows, ignore_index=True)
        outp = REG_DIR / "edges_long.parquet"
        edges_long.to_parquet(outp, index=False)
        print(f"↳ edges_long written → {outp}  (warning: large file)")

    print(f"→ valid subjects for regression: {row_valid.sum()}/{len(row_valid)}")


# ------------------------ Entrypoint ------------------------
if __name__ == "__main__":
    main()
