# 4_statistics/2_run_ttest.py
"""
Run Welch two-sample t-tests (ASD vs HC) on network scores and apply FDR
correction within families (intra / inter). Outputs a CSV where HC/ASD are
shown as signed 'mean ± sd' (4 decimals), plus dispersion stats and
adjusted p-values.

Inputs (in stage_artifacts/4_statistics):
- network_scores.csv  (from 1_build_network_scores.py)

Outputs (in stage_artifacts/4_statistics):
- network_tests.csv

Notes
-----
- Group order for t-tests is ASD_vs_HC: t > 0 means ASD > HC.
- HC/ASD columns preserve the sign (no absolute or '-0.0000' normalization).
"""

# ------------------------ Imports ------------------------
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

from config import ALPHA, FDR_METHOD, stage_artifacts

# ------------------------ Artifacts ------------------------
ART_STAT = stage_artifacts("4_statistics")

# ------------------------ Formatting helpers ------------------------
def fmt_signed(v: float, nd: int = 4) -> str:
    """Format a float with nd decimals, preserving the sign."""
    return f"{v:.{nd}f}"

def fmt_p_disp(p: float) -> str:
    """Pretty p display: '<0.001' or 3 decimals."""
    return "<0.001" if p < 1e-3 else f"{p:.3f}"

# ------------------------ Core ------------------------
def run_for(scores_csv: Path, out_csv: Path, nd: int) -> None:
    """
    Run Welch t-tests (ASD vs HC) per (type, pair) and apply FDR within 'type'.
    Saves HC/ASD as signed 'mean ± sd' strings with nd decimals.
    """
    df = pd.read_csv(scores_csv)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")

    rows = []
    for (typ, pair), g in df.groupby(["type", "pair"], sort=False):
        a = g.loc[g["group"] == "HC", "score"].to_numpy(dtype=float)   # HC
        b = g.loc[g["group"] == "ASD", "score"].to_numpy(dtype=float)  # ASD
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        if len(a) < 3 or len(b) < 3:
            continue

        # Welch t-test (ASD vs HC) -> t > 0 means ASD > HC
        t, p = ttest_ind(b, a, equal_var=False)

        # Summary stats
        m_hc  = float(np.nanmean(a)); sd_hc  = float(np.nanstd(a, ddof=1)); var_hc  = float(sd_hc**2)
        m_asd = float(np.nanmean(b)); sd_asd = float(np.nanstd(b, ddof=1)); var_asd = float(sd_asd**2)

        # Display strings: signed mean ± sd
        hc_disp  = f"{fmt_signed(m_hc, nd)} ± {fmt_signed(sd_hc, nd)}"
        asd_disp = f"{fmt_signed(m_asd, nd)} ± {fmt_signed(sd_asd, nd)}"

        rows.append({
            "type": typ,
            "Connection": pair,
            "HC": hc_disp,
            "ASD": asd_disp,
            # stats inférentielles
            "t value": float(t),
            "p_raw": float(p),
            "group_order_for_ttest": "ASD_vs_HC",
        })

    out = pd.DataFrame(rows)
    if out.empty:
        print(f"[!] No valid tests run for {scores_csv}")
        return

    # FDR within each family (intra / inter)
    out["p_FDR"] = np.nan
    for typ, g in out.groupby("type"):
        _, p_adj, *_ = multipletests(g["p_raw"].values, alpha=ALPHA, method=FDR_METHOD)
        out.loc[g.index, "p_FDR"] = p_adj

    # Pretty p-value
    out["p-value"] = out["p_FDR"].apply(fmt_p_disp)

    # Logs
    n_total = len(out)
    n_sig_uncorr = int((out["p_raw"] < ALPHA).sum())
    n_sig_fdr = int((out["p_FDR"] < ALPHA).sum())
    n_pos = int((out["t value"] > 0).sum())
    n_neg = int((out["t value"] < 0).sum())

    # Persist (drop raw p; keep adjusted + pretty)
    out = out.sort_values(["type", "Connection"]).drop(columns=["p_raw"])
    out.to_csv(out_csv, index=False)

    print(f"\n=== Results for {scores_csv.name} ===")
    print(f"Total tested blocks : {n_total}")
    print(f"  - Significant before FDR : {n_sig_uncorr}")
    print(f"  - Significant after FDR  : {n_sig_fdr}")
    print(f"  - t>0 : {n_pos}   t<0 : {n_neg}")
    print(f"✓ Saved → {out_csv}")

# ------------------------ Entrypoint ------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--decimals",
        type=int,
        default=4,
        help="Decimals for 'mean ± sd' display. Default: 4.",
    )
    args = parser.parse_args()

    s = ART_STAT / "network_scores.csv"
    if s.exists():
        run_for(s, ART_STAT / "network_tests.csv", nd=args.decimals)
    else:
        print(f"[!] Missing input: {s}")

if __name__ == "__main__":
    main()
