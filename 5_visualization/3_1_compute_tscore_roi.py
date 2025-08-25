# 4_statistics/1_compute_tscore_roi.py

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

from config import stage_artifacts, ALPHA, FDR_METHOD, USE_EQUAL_VAR

# ------------------------ Artifact folders ------------------------
ART_STAT = stage_artifacts("4_statistics")
ART_STAT.mkdir(parents=True, exist_ok=True)

# NOTE: fixed the path join here (no string division).
ART_VIZ = stage_artifacts("5_visualization")
ART_SURF = ART_VIZ / "surfaces_roi_t"
ART_SURF.mkdir(parents=True, exist_ok=True)

# ------------------------ Inputs ------------------------
CONN_BASE   = stage_artifacts("2_connectivity_extraction") / "regressed"
MEMMAP_PATH = CONN_BASE / "conn_resid.memmap"
LABELS_PATH = CONN_BASE / "labels_valid.npy"
NROIS_PATH  = CONN_BASE / "n_rois.txt"
TRIU_I_PATH = CONN_BASE / "triu_i.npy"
TRIU_J_PATH = CONN_BASE / "triu_j.npy"
X_VEC_PATH  = CONN_BASE / "connectivity_residuals.npy"  # (n_subj, n_edges) – optional
ROI_HO_CSV = stage_artifacts("3_annotation_networks") / "roi_to_HO.csv"

ROI_YEO_CSV = stage_artifacts("3_annotation_networks") / "roi_to_yeo.csv"

# --- Helpers for safe label handling (put near the top) ---
def _norm_label(x) -> str:
    """Return a clean string label, replacing NaN/None/empty by 'Unknown'."""
    if x is None:
        return "Unknown"
    try:
        # pandas NaN/NA check
        import pandas as pd
        if pd.isna(x):
            return "Unknown"
    except Exception:
        pass
    s = str(x).strip()
    return s if s else "Unknown"

def _pair_label(a, b) -> str:
    """Create a symmetric 'A ↔ B' pair with safe normalization."""
    A, B = _norm_label(a), _norm_label(b)
    return " ↔ ".join(sorted([A, B]))


# ------------------------ Helpers ------------------------
def load_vectorized_or_memmap(n_subj: int, n_rois: int, iu):
    """
    Load edge matrix X (subjects x edges). Prefer pre-vectorized .npy if present,
    otherwise vectorize the memmap on the fly.

    Parameters
    ----------
    n_subj : int
    n_rois : int
    iu : tuple(np.ndarray, np.ndarray)
        Upper-triangle indices (k=1).

    Returns
    -------
    np.ndarray
        X of shape (n_subj, n_edges).
    """
    if X_VEC_PATH.exists():
        X = np.load(X_VEC_PATH)
        assert X.shape == (n_subj, iu[0].size), f"X shape {X.shape} != {(n_subj, iu[0].size)}"
        return X

    assert MEMMAP_PATH.exists(), f"Not found: {MEMMAP_PATH}"
    conn = np.memmap(MEMMAP_PATH, dtype="float32", mode="r", shape=(n_subj, n_rois, n_rois))
    X = np.empty((n_subj, iu[0].size), dtype="float32")
    for s in range(n_subj):
        X[s] = conn[s][iu]
    return X


# ------------------------ Main ------------------------
def main():
    # --- Load labels and n_rois ---
    assert LABELS_PATH.exists(), f"Not found: {LABELS_PATH}"
    labels = np.load(LABELS_PATH)  # 1=ASD, 2=HC
    n_subj = int(labels.shape[0])

    assert NROIS_PATH.exists(), f"Not found: {NROIS_PATH}"
    n_rois = int(Path(NROIS_PATH).read_text().strip())

    # upper-triangle indices
    if TRIU_I_PATH.exists() and TRIU_J_PATH.exists():
        iu = (np.load(TRIU_I_PATH), np.load(TRIU_J_PATH))
    else:
        iu = np.triu_indices(n_rois, k=1)

    # --- Edge matrix X ---
    X = load_vectorized_or_memmap(n_subj, n_rois, iu)

    # --- Groups ---
    X_asd = X[labels == 1]
    X_hc  = X[labels == 2]

    # --- Welch (or Student) t-tests per edge ---
    t_vals, p_unc = ttest_ind(
        X_asd, X_hc, axis=0,
        equal_var=bool(USE_EQUAL_VAR),
        nan_policy="omit"
    )

    # --- FDR correction across all edges ---
    rejected, p_fdr, _, _ = multipletests(p_unc, alpha=ALPHA, method=FDR_METHOD)

    # --- Build symmetric ROI×ROI maps ---
    t_map    = np.zeros((n_rois, n_rois), dtype=np.float32)
    pmap_unc = np.ones((n_rois, n_rois), dtype=np.float32)
    pmap_fdr = np.ones((n_rois, n_rois), dtype=np.float32)
    sig_mask = np.zeros((n_rois, n_rois), dtype=bool)

    t_map[iu]    = t_vals
    pmap_unc[iu] = p_unc
    pmap_fdr[iu] = p_fdr
    sig_mask[iu] = rejected

    # symmetrize
    t_map     = t_map + t_map.T
    pmap_unc  = pmap_unc + pmap_unc.T - np.eye(n_rois, dtype=float)
    pmap_fdr  = pmap_fdr + pmap_fdr.T - np.eye(n_rois, dtype=float)
    sig_mask  = sig_mask | sig_mask.T

    # --- Save main arrays ---
    out_dir = ART_SURF / "edge_ttests_fdr"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / "t_map.npy", t_map)
    np.save(out_dir / "p_unc_map.npy", pmap_unc)
    np.save(out_dir / "p_fdr_map.npy", pmap_fdr)
    np.save(out_dir / "sig_mask.npy", sig_mask)
    np.save(out_dir / "triu_i.npy", iu[0])
    np.save(out_dir / "triu_j.npy", iu[1])
    (out_dir / "n_rois.txt").write_text(str(n_rois))

    # --- Significant edges table (no averaging) ---
    i_e, j_e = iu
    df = pd.DataFrame({
        "roi_i": (i_e + 1).astype(np.int32),
        "roi_j": (j_e + 1).astype(np.int32),
        "tval":  t_vals.astype(np.float32),
        "p_unc": p_unc.astype(np.float32),
        "p_fdr": p_fdr.astype(np.float32),
        "significant": rejected.astype(bool),
    })
    df_sig = df[df["significant"]].copy()

    # --- Optional enrichment with Yeo7/17 labels ---
    if ROI_YEO_CSV.exists():
        m = pd.read_csv(ROI_YEO_CSV)
        m = m.rename(columns={
            "ROI_number": "roi",
            "Yeo7_id": "yeo7_id",   "Yeo7_name": "yeo7_name",
            "Yeo17_id": "yeo17_id", "Yeo17_name": "yeo17_name",
        })

        df_sig = df_sig.merge(
            m.rename(columns={
                "roi": "roi_i",
                "yeo7_id": "yeo7_i", "yeo7_name": "yeo7_i_name",
                "yeo17_id": "yeo17_i", "yeo17_name": "yeo17_i_name",
            }),
            on="roi_i", how="left"
        )
        df_sig = df_sig.merge(
            m.rename(columns={
                "roi": "roi_j",
                "yeo7_id": "yeo7_j", "yeo7_name": "yeo7_j_name",
                "yeo17_id": "yeo17_j", "yeo17_name": "yeo17_j_name",
            }),
            on="roi_j", how="left"
        )

    if ROI_HO_CSV.exists():
        # force dtype to avoid numeric NaN sneaking in as floats
        ho = pd.read_csv(ROI_HO_CSV, dtype={"ROI_number": int, "HO_label": "string"}) \
            .rename(columns={"ROI_number": "roi", "HO_label": "HO_name"})

        df_sig = df_sig.merge(
            ho.rename(columns={"roi": "roi_i", "HO_name": "HO_i_name"}),
            on="roi_i", how="left"
        )
        df_sig = df_sig.merge(
            ho.rename(columns={"roi": "roi_j", "HO_name": "HO_j_name"}),
            on="roi_j", how="left"
        )

        # Safe pair (handles NaN/None/empty)
        df_sig["HO_pair"] = [
            _pair_label(a, b) for a, b in zip(df_sig.get("HO_i_name"), df_sig.get("HO_j_name"))
        ]



        def pair(a, b):
            a, b = str(a or "Unknown"), str(b or "Unknown")
            return " ↔ ".join(sorted([a, b]))

        df_sig["yeo7_pair"]  = [pair(a, b) for a, b in zip(df_sig["yeo7_i_name"],  df_sig["yeo7_j_name"])]
        df_sig["yeo7_kind"]  = np.where(df_sig["yeo7_i_name"] == df_sig["yeo7_j_name"], "intra", "inter")
        df_sig["yeo17_pair"] = [pair(a, b) for a, b in zip(df_sig["yeo17_i_name"], df_sig["yeo17_j_name"])]
        df_sig["yeo17_kind"] = np.where(df_sig["yeo17_i_name"] == df_sig["yeo17_j_name"], "intra", "inter")

    df_sig = df_sig.sort_values("tval", ascending=False)

    # --- Write edge tables ---
    df_sig.to_parquet(out_dir / "edges_significant.parquet", index=False)
    df_sig.to_csv(out_dir / "edges_significant.csv", index=False)

    # --- Global summary (counts only; no means) ---
    n_edges_tested = int(t_vals.size)
    n_edges_sig = int(df_sig.shape[0])
    n_asd = int((labels == 1).sum())
    n_hc  = int((labels == 2).sum())
    pd.DataFrame([{
        "alpha": float(ALPHA),
        "fdr_method": str(FDR_METHOD),
        "use_equal_var": bool(USE_EQUAL_VAR),
        "n_subjects_total": int(labels.size),
        "n_subjects_ASD": n_asd,
        "n_subjects_HC": n_hc,
        "n_rois": int(n_rois),
        "n_edges_tested": n_edges_tested,
        "n_edges_significant": n_edges_sig,
        "prop_significant": float(n_edges_sig / max(1, n_edges_tested)),
    }]).to_csv(out_dir / "summary_global.csv", index=False)

    # --- Summaries (counts only) ---
    def summarize_intra(df_edges: pd.DataFrame, col_name: str) -> pd.DataFrame:
        if df_edges.empty:
            return pd.DataFrame(columns=["network", "n_edges", "n_unique_rois"])
        rows = []
        for net, grp in df_edges.groupby(col_name):
            rois = np.unique(np.concatenate([grp.roi_i.values, grp.roi_j.values]))
            rows.append({"network": net, "n_edges": int(len(grp)), "n_unique_rois": int(len(rois))})
        return pd.DataFrame(rows).sort_values(["n_edges", "network"], ascending=[False, True])

    def summarize_inter(df_edges: pd.DataFrame, pair_col: str) -> pd.DataFrame:
        if df_edges.empty:
            return pd.DataFrame(columns=["pair", "n_edges", "n_unique_rois"])
        rows = []
        for p, grp in df_edges.groupby(pair_col):
            rois = np.unique(np.concatenate([grp.roi_i.values, grp.roi_j.values]))
            rows.append({"pair": p, "n_edges": int(len(grp)), "n_unique_rois": int(len(rois))})
        return pd.DataFrame(rows).sort_values(["n_edges", "pair"], ascending=[False, True])

    # ROI degree (over all significant edges)
    def roi_degree(df_edges: pd.DataFrame, n_rois: int, map_df: pd.DataFrame | None = None) -> pd.DataFrame:
        deg = np.zeros(n_rois, dtype=np.int32)
        if not df_edges.empty:
            i = df_edges["roi_i"].to_numpy(int) - 1
            j = df_edges["roi_j"].to_numpy(int) - 1
            for a, b in zip(i, j):
                if 0 <= a < n_rois:
                    deg[a] += 1
                if 0 <= b < n_rois:
                    deg[b] += 1
        out = pd.DataFrame({"ROI_number": np.arange(1, n_rois + 1, dtype=int), "deg_sig": deg})
        if map_df is not None:
            out = out.merge(
                map_df[["ROI_number", "Yeo7_id", "Yeo7_name", "Yeo17_id", "Yeo17_name"]],
                on="ROI_number", how="left"
            )
        return out.sort_values("deg_sig", ascending=False)

    # Write Yeo7 / Yeo17 summaries if mapping columns are present
    if {"yeo7_kind", "yeo7_i_name", "yeo7_pair"}.issubset(df_sig.columns):
        df7_intra = df_sig[df_sig["yeo7_kind"] == "intra"]
        df7_inter = df_sig[df_sig["yeo7_kind"] == "inter"]
        summarize_intra(df7_intra, "yeo7_i_name").to_csv(out_dir / "summary_yeo7_intra.csv", index=False)
        summarize_inter(df7_inter, "yeo7_pair").to_csv(out_dir / "summary_yeo7_inter.csv", index=False)

    if {"yeo17_kind", "yeo17_i_name", "yeo17_pair"}.issubset(df_sig.columns):
        df17_intra = df_sig[df_sig["yeo17_kind"] == "intra"]
        df17_inter = df_sig[df_sig["yeo17_kind"] == "inter"]
        summarize_intra(df17_intra, "yeo17_i_name").to_csv(out_dir / "summary_yeo17_intra.csv", index=False)
        summarize_inter(df17_inter, "yeo17_pair").to_csv(out_dir / "summary_yeo17_inter.csv", index=False)

    # ROI degree CSV (with Yeo mapping if available)
    map_df = pd.read_csv(ROI_YEO_CSV) if ROI_YEO_CSV.exists() else None
    roi_degree(df_sig, n_rois, map_df).to_csv(out_dir / "roi_degree_all.csv", index=False)

    print(
        f"✓ Files written → {out_dir} : "
        f"edges_significant.csv, summary_global.csv, summary_yeo7_*.csv, "
        f"summary_yeo17_*.csv, roi_degree_all.csv, and NumPy maps."
    )


# ------------------------ Entrypoint ------------------------
if __name__ == "__main__":
    main()
