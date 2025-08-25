# 4_statistics/1_build_network_scores.py

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

# ------------------------ Matplotlib backend ------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind
from config import stage_artifacts

# ------------------------ Artifacts ------------------------
ART_CONN = stage_artifacts("2_connectivity_extraction")
REG_DIR  = ART_CONN / "regressed"
ART_ANN  = stage_artifacts("3_annotation_networks")
ART_OUT  = stage_artifacts("4_statistics")
ART_OUT.mkdir(parents=True, exist_ok=True)
FIG_DIR  = ART_OUT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------ Abbreviations ------------------------
ABBR = {
    "Default Mode": "DMN",
    "Dorsal Attention": "DAN",
    "Frontoparietal Control": "FPN",
    "Limbic": "LIM",
    "Salience / Ventral Attention": "SAL",
    "Somatomotor": "SMN",
    "Visual": "VIS",
}

def short(x: str) -> str:
    return ABBR.get(x, x)

# ------------------------ Loading ------------------------

def _open_residual_memmap():
    conn_mm = REG_DIR / "conn_resid.memmap"
    assert conn_mm.exists(), f"Not found: {conn_mm}"
    n_rois = int((REG_DIR / "n_rois.txt").read_text().strip())
    subjects = np.load(REG_DIR / "subjects_valid.npy", allow_pickle=True)
    n_subj = len(subjects)
    conn = np.memmap(conn_mm, dtype="float32", mode="r", shape=(n_subj, n_rois, n_rois))
    return conn, subjects, n_rois


def _load_groups(subjects):
    labels = np.load(REG_DIR / "labels_valid.npy")
    # ABIDE: 1=ASD, 2=HC
    grp = np.where(labels == 1, "ASD", "HC")
    return np.asarray(subjects, dtype=object), grp


def _load_yeo_mapping(n_rois: int, yeo: str = "7"):
    """
    Load ROI -> Yeo network mapping (7 or 17).
    Expected: roi_to_yeo.csv with columns:
      ROI_number, Yeo7_id, Yeo7_name, Yeo17_id, Yeo17_name
    Returns
      - nets : sorted list of network names (without "Unknown")
      - roi_idx : dict {network_name: numpy_roi_indices}
    """
    assert yeo in {"7", "17"}, "yeo must be '7' or '17'"
    col = "Yeo7_name" if yeo == "7" else "Yeo17_name"
    y = pd.read_csv(ART_ANN / "roi_to_yeo.csv")
    y = y[["ROI_number", col]].rename(columns={"ROI_number": "roi", col: "net"})
    y["net"] = y["net"].fillna("Unknown")
    assert len(y) >= n_rois, "roi_to_yeo.csv does not cover all ROIs"
    nets = [n for n in y["net"].unique().tolist() if n != "Unknown"]
    roi_idx = {net: y.loc[y["net"] == net, "roi"].to_numpy(dtype=int) - 1 for net in nets}
    return sorted(nets), roi_idx

# ------------------------ Sparsity masks ------------------------

def build_group_mask(conn: np.ndarray, sparsity: float, keep: str = "abs") -> np.ndarray:
    """Global mask (same edges kept for all)."""
    n_subj, n_rois, _ = conn.shape
    iu = np.triu_indices(n_rois, 1)
    vals = conn[:, iu[0], iu[1]]  # (n_subj, n_edges)

    if keep == "pos":
        magn = np.mean(vals * (vals > 0), axis=0)
    elif keep == "neg":
        magn = np.mean(-vals * (vals < 0), axis=0)
    else:  # abs
        magn = np.mean(np.abs(vals), axis=0)

    k = max(1, int(np.floor(sparsity * len(magn))))
    thr = np.partition(magn, -k)[-k]
    keep_edges = magn >= thr

    mask = np.zeros((n_rois, n_rois), dtype=bool)
    mask[iu] = keep_edges
    mask |= mask.T
    return mask


def build_subject_mask(M: np.ndarray, sparsity: float, keep: str = "abs") -> np.ndarray:
    """Subject-specific mask (top-k% on |z|, pos, or neg)."""
    n_rois = M.shape[0]
    iu = np.triu_indices(n_rois, 1)
    v = M[iu]

    if keep == "pos":
        sel = v > 0
        magn = v.copy()
    elif keep == "neg":
        sel = v < 0
        magn = -v.copy()
    else:
        sel = np.ones_like(v, dtype=bool)
        magn = np.abs(v)

    n_edges = int(sel.sum())
    if n_edges == 0:
        return np.zeros_like(M, dtype=bool)

    k = max(1, int(np.floor(sparsity * n_edges)))
    thr = np.partition(magn[sel], -k)[-k]
    keep_upper = sel & (magn >= thr)

    mask = np.zeros((n_rois, n_rois), dtype=bool)
    mask[iu] = keep_upper
    mask |= mask.T
    np.fill_diagonal(mask, False)
    return mask

# ------------------------ Mean scores ------------------------

def build_scores_mean(conn, subjects, grp, nets, roi_idx, sparsity, keep, mask_mode):
    rows = []
    n_subj, n_rois, _ = conn.shape
    nets_sorted = sorted(nets)
    global_mask = None

    if sparsity > 0 and mask_mode == "global":
        global_mask = build_group_mask(conn, sparsity, keep)
        n_total = n_rois * (n_rois - 1) // 2
        n_kept = int(np.triu(global_mask, 1).sum())
        print(f"→ Kept edges (global): {n_kept}/{n_total} ({100*n_kept/n_total:.1f}%)")

    for s in range(n_subj):
        M = conn[s].copy()
        np.fill_diagonal(M, np.nan)

        if sparsity > 0:
            mask_s = global_mask if mask_mode == "global" else build_subject_mask(M, sparsity, keep)
            if mask_mode == "subject":
                n_total = n_rois * (n_rois - 1) // 2
                n_kept = int(np.triu(mask_s, 1).sum())
                if s == 0:
                    print(f"→ Kept edges (subject): {100*n_kept/n_total:.1f}% per subject")
            M[~mask_s] = np.nan

        # Intra-network
        for net in nets_sorted:
            idx = roi_idx[net]
            if idx.size < 2:
                continue
            sub = M[np.ix_(idx, idx)]
            vals = sub[np.triu_indices_from(sub, 1)]
            score = float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan
            rows.append((subjects[s], grp[s], "intra", short(net), score))

        # Inter-network
        for i, a in enumerate(nets_sorted):
            for b in nets_sorted[i+1:]:
                ia, ib = roi_idx[a], roi_idx[b]
                sub = M[np.ix_(ia, ib)]
                score = float(np.nanmean(sub)) if np.isfinite(sub).any() else np.nan
                rows.append((subjects[s], grp[s], "inter", f"{short(a)}-{short(b)}", score))

    return pd.DataFrame(rows, columns=["subject","group","type","pair","score"])

# ------------------------ PCA utilities ------------------------

def pca_from_edges_matrix(X: np.ndarray):
    """
    PCA via SVD on X (n_subj x n_edges) with column centering.
    Returns:
      scores : (n_subj, n_comp) = U * S
      var_exp: (n_comp,) proportions of explained variance
      Vt     : (n_comp, n_edges) loadings
    """
    col_mean = np.nanmean(X, axis=0, keepdims=True)
    Xc = X - col_mean
    Xc = np.where(np.isfinite(Xc), Xc, 0.0)

    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    s2 = S**2
    var_exp = s2 / s2.sum() if s2.sum() > 0 else np.zeros_like(s2)
    scores = U * S
    return scores, var_exp, Vt


def _zscore_cols(A: np.ndarray, ddof: int = 1) -> np.ndarray:
    m = A.mean(axis=0, keepdims=True)
    s = A.std(axis=0, ddof=ddof, keepdims=True)
    return (A - m) / (s + 1e-12)


def _align_scores_to_anchor(scores: np.ndarray, anchor: np.ndarray, k_align: int) -> None:
    """Align sign of the first k components to an anchor vector (in-place)."""
    if not np.isfinite(anchor).any() or np.std(anchor) == 0:
        return
    for j in range(min(k_align, scores.shape[1])):
        r = np.corrcoef(scores[:, j], anchor)[0, 1]
        if np.isfinite(r) and r < 0:
            scores[:, j] *= -1.0

# ------------------------ PCA scores ------------------------

def build_scores_pca(conn, subjects, grp, nets, roi_idx, sparsity, keep, mask_mode,
                     n_pca: int = 1, report_top=5, make_report=True,
                     anchor_mode: str = "signed"):
    """
    Builds:
      - network_scores.csv : subject score per block (PC1 or Composite_PC1..PCk)
      - (optional) network_pca_report.csv : varExp PC1..PCn + t-tests PC1..PCn + composite t-test
      - Figures: explained variance / violin plots / cumulative explained variance
    Sign anchoring (anchor_mode):
      - 'signed' (default): higher = more positive mean connectivity (signed mean of block)
      - 'abs'              : higher = stronger connectivity magnitude (mean of |edges|)
      - 'group'            : higher = higher values in ASD (ASD=+1, HC=-1)
    """
    n_subj, n_rois, _ = conn.shape
    nets_sorted = sorted(nets)
    global_mask = None
    if sparsity > 0 and mask_mode == "global":
        global_mask = build_group_mask(conn, sparsity, keep)
        n_total = n_rois * (n_rois - 1) // 2
        n_kept = int(np.triu(global_mask, 1).sum())
        print(f"→ Kept edges (global): {n_kept}/{n_total} ({100*n_kept/n_total:.1f}%)")

    rows_scores, rows_report = [], []

    # Blocks
    blocks = [("intra", short(n)) for n in nets_sorted]
    for i, a in enumerate(nets_sorted):
        for b in nets_sorted[i+1:]:
            blocks.append(("inter", f"{short(a)}-{short(b)}"))

    for typ, pair in blocks:
        rows_X = []
        for s in range(n_subj):
            M = conn[s].copy()
            np.fill_diagonal(M, np.nan)
            if sparsity > 0:
                mask_s = global_mask if mask_mode == "global" else build_subject_mask(M, sparsity, keep)
                if mask_mode == "subject" and s == 0:
                    n_total = n_rois * (n_rois - 1) // 2
                    n_kept = int(np.triu(mask_s, 1).sum())
                    print(f"→ Kept edges (subject): {100*n_kept/n_total:.1f}% per subject")
                M[~mask_s] = np.nan

            if typ == "intra":
                net = pair
                long_net = [k for k, v in ABBR.items() if v == net]
                long_net = long_net[0] if long_net else net
                idx = roi_idx[long_net]
                sub = M[np.ix_(idx, idx)]
                v = sub[np.triu_indices_from(sub, 1)]
            else:
                a, b = pair.split("-")
                long_a = [k for k, v in ABBR.items() if v == a]; long_a = long_a[0] if long_a else a
                long_b = [k for k, v in ABBR.items() if v == b]; long_b = long_b[0] if long_b else b
                ia, ib = roi_idx[long_a], roi_idx[long_b]
                v = M[np.ix_(ia, ib)].reshape(-1)

            rows_X.append(v)

        X = np.vstack(rows_X)
        col_mean = np.nanmean(X, axis=0)
        inds = np.where(~np.isfinite(X))
        X[inds] = np.take(col_mean, inds[1])

        if anchor_mode == "group":
            anchor = np.where(grp == "ASD", 1.0, -1.0)
        elif anchor_mode == "abs":
            anchor = np.mean(np.abs(X), axis=1)
        else:
            anchor = np.mean(X, axis=1)

        scores, var_exp, _ = pca_from_edges_matrix(X)
        if scores.size == 0:
            for sid, g in zip(subjects, grp):
                rows_scores.append((sid, g, typ, pair, np.nan, "NA", anchor_mode))
            continue

        max_k = scores.shape[1]
        k_sel = int(max(1, min(n_pca, max_k)))
        _align_scores_to_anchor(scores, anchor, k_sel)

        PCs = scores[:, :k_sel]
        PCs_z = _zscore_cols(PCs)

        w = var_exp[:k_sel].copy()
        w = w / (float(w.sum()) + 1e-12) if w.sum() > 0 else np.ones_like(w)/len(w)

        comp = (PCs_z @ w.reshape(-1, 1)).ravel()
        comp_z = (comp - comp.mean()) / (comp.std(ddof=1) + 1e-12)

        pc1_z = (scores[:, 0] - scores[:, 0].mean()) / (scores[:, 0].std(ddof=1) + 1e-12)

        score_kind = "PC1" if k_sel == 1 else f"Composite_PC1..PC{k_sel}"
        final_score = pc1_z if k_sel == 1 else comp_z
        for sid, g, sc in zip(subjects, grp, final_score):
            rows_scores.append((sid, g, typ, pair, float(sc), score_kind, anchor_mode))

        if make_report:
            report = {"type": typ, "pair": pair, "n_edges": int(X.shape[1])}
            for j in range(k_sel):
                report[f"var_PC{j+1}"] = float(var_exp[j])
            for j in range(k_sel):
                pc_z = (scores[:, j] - scores[:, j].mean()) / (scores[:, j].std(ddof=1) + 1e-12)
                a = pc_z[grp == "HC"]; b = pc_z[grp == "ASD"]
                t, p = ttest_ind(b, a, equal_var=False)
                report[f"t_PC{j+1}"] = float(t); report[f"p_PC{j+1}"] = float(p)
            a = final_score[grp == "HC"]; b = final_score[grp == "ASD"]
            t, p = ttest_ind(b, a, equal_var=False)
            report["t_Final"] = float(t); report["p_Final"] = float(p)
            rows_report.append(report)

    # === Write CSVs ===
    df_scores = pd.DataFrame(rows_scores, columns=["subject","group","type","pair","score","score_kind","anchor_mode"])
    out_scores = ART_OUT / "network_scores.csv"
    df_scores.to_csv(out_scores, index=False)
    print(f"✓ Wrote network scores → {out_scores}  [rows={len(df_scores)}]  ({df_scores['score_kind'].unique().tolist()})")

    if make_report and len(rows_report) > 0:
        df_rep = pd.DataFrame(rows_report)
        out_rep = ART_OUT / "network_pca_report.csv"
        df_rep.to_csv(out_rep, index=False)
        print(f"✓ Wrote PCA report → {out_rep}")

        # === Figures ===
        # 1) Barplot: PC1 explained variance by block
        df_rep_plot = df_rep.copy()
        if "var_PC1" not in df_rep_plot.columns:
            df_rep_plot["var_PC1"] = 0.0
        df_rep_plot["block"] = df_rep_plot["type"].str.upper() + " • " + df_rep_plot["pair"]
        df_rep_plot = df_rep_plot.sort_values("var_PC1", ascending=False)

        fig = plt.figure(figsize=(10, 6))
        plt.barh(df_rep_plot["block"], df_rep_plot["var_PC1"])
        plt.gca().invert_yaxis()
        plt.xlabel("Explained variance by PC1")
        plt.title("PCA by block — PC1 explained variance")
        plt.tight_layout()
        p1 = FIG_DIR / "pca_var_explained_PC1_by_blocks.png"
        fig.savefig(str(p1), dpi=200, bbox_inches="tight"); plt.close(fig)
        print(f"✓ Figure → {p1}")

        # 2) Stacked bars: cumulative explained variance (PC1..PCk)
        pc_cols = sorted([int(c.replace("var_PC","")) for c in df_rep.columns if c.startswith("var_PC")])
        if len(pc_cols) == 0:
            pc_cols = [1]
        k_plot = int(min(n_pca, max(pc_cols)))
        pcs = [f"var_PC{j}" for j in range(1, k_plot+1)]
        for c in pcs:
            if c not in df_rep.columns:
                df_rep[c] = 0.0

        df_cum = df_rep[["type", "pair"] + pcs].copy()
        df_cum[pcs] = df_cum[pcs].fillna(0.0)
        df_cum["cumvar_1..k"] = df_cum[pcs].sum(axis=1)
        df_cum["block"] = df_cum["type"].str.upper() + " • " + df_cum["pair"]
        df_cum = df_cum.sort_values("cumvar_1..k", ascending=False).head(report_top)
        labels = df_cum["block"].tolist()

        fig = plt.figure(figsize=(10, max(4, 0.35 * len(labels) + 2)))
        left = np.zeros(len(df_cum), dtype=float)
        for j, pc in enumerate(pcs, start=1):
            vals = df_cum[pc].values
            plt.barh(labels, vals, left=left, label=f"PC{j}")
            left += vals
        for y, v in enumerate(left):
            plt.text(v, y, f"{v:.2f}", va="center", ha="left", fontsize=8)
        plt.gca().invert_yaxis()
        plt.xlabel(f"Cumulative explained variance (PC1..PC{k_plot})")
        plt.title(f"PCA by block — cumulative (PC1..PC{k_plot})")
        plt.legend(title="Components", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        p3 = FIG_DIR / f"pca_cumvar_PC1..PC{k_plot}_top_blocks.png"
        fig.savefig(str(p3), dpi=200, bbox_inches="tight"); plt.close(fig)
        print(f"✓ Figure → {p3}")

    return df_scores

# ------------------------ Main ------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sparsity", type=float, default=0.2,
                        help="Sparsity proportion (0.2 = 20%).")
    parser.add_argument("--keep", choices=["abs","pos","neg"], default="abs",
                        help="Edge ranking criterion for sparsity.")
    parser.add_argument("--mask-mode", choices=["global","subject","none"], default="global",
                        help="Global sparsity (same mask), per-subject, or none.")
    parser.add_argument("--score-mode", choices=["mean","pca"], default="pca",
                        help="Block summary: simple mean or PCA (PC1/composite).")
    parser.add_argument("--pca-report", action="store_true",
                        help="Write network_pca_report.csv and figures.")
    parser.add_argument("--report-top", type=int, default=5,
                        help="Number of top blocks for violins and cumulative plots.")
    parser.add_argument("--n-pca", type=int, default=3,
                        help="Number of PCA components for the score (1=PC1; >=2=weighted composite).")
    parser.add_argument("--anchor-mode", choices=["signed","abs","group"], default="signed",
                        help="Sign anchoring: 'signed'=signed mean (higher=more positive), 'abs'=magnitude, 'group'=ASD>HC.")
    parser.add_argument("--yeo", choices=["7","17"], default="7",
                        help="Yeo parcellation to annotate ROIs (7 or 17).")
    args = parser.parse_args()

    conn, subjects, n_rois = _open_residual_memmap()
    subjects, grp = _load_groups(subjects)
    nets, roi_idx = _load_yeo_mapping(n_rois, yeo=args.yeo)

    print(f"→ Connectivity: {len(subjects)} subjects, {n_rois} ROIs")
    print(f"→ Annotation: Yeo{args.yeo}")
    print(f"→ Sparsity: {args.sparsity:.2f} ({args.keep}) [mode={args.mask_mode}]")
    if args.score_mode == "pca":
        print(f"→ Score mode: {args.score_mode} | n_pca={args.n_pca} | anchor={args.anchor_mode}")
    else:
        print(f"→ Score mode: {args.score_mode}")

    sparsity = 0.0 if args.mask_mode == "none" else float(args.sparsity)

    if args.score_mode == "mean":
        df = build_scores_mean(conn, subjects, grp, nets, roi_idx,
                               sparsity, args.keep, args.mask_mode)
        out = ART_OUT / "network_scores.csv"
        df.to_csv(out, index=False)
        print(f"✓ Wrote mean scores → {out}  [rows={len(df)}]")
    else:
        build_scores_pca(conn, subjects, grp, nets, roi_idx,
                         sparsity, args.keep, args.mask_mode,
                         n_pca=args.n_pca, report_top=args.report_top,
                         make_report=args.pca_report, anchor_mode=args.anchor_mode)


if __name__ == "__main__":
    main()
