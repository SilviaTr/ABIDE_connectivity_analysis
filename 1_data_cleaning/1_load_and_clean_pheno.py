# 1_data_cleaning/1_load_and_clean_pheno.py

import pandas as pd
import numpy as np
from pathlib import Path
from config import PHENO_FILE, stage_artifacts, DATA_DIR

# ------------------------ QC Parameters ------------------------
N_MIN_CELL = 10          # warn if cells < N_MIN_CELL in SITE×DX
DELTA_PROP_SEX = 0.15    # warn if difference in male proportion >= 0.15
ALPHA = 0.05             # p-value threshold (kept for reference)
AGE_BINS = [0, 12, 18, 30, 50, 120]             # age bins
AGE_LABELS = ["≤12", "13–18", "19–30", "31–50", "≥51"]
# ---------------------------------------------------------------

STAGE = Path(__file__).resolve().parent.name
ART = stage_artifacts(STAGE)
ART.mkdir(parents=True, exist_ok=True)


# ------------------------ Utilities ------------------------
def try_import_scipy():
    """
    Try to import scipy.stats. Return the module if available, else None.
    This keeps the script robust on environments without SciPy.
    """
    try:
        from scipy import stats
        return stats
    except Exception:
        return None


def welch_ttest(x, y):
    """
    Welch's t-test (manual implementation without SciPy).

    Parameters
    ----------
    x, y : array-like
        Samples.

    Returns
    -------
    t : float
        t statistic.
    df : float
        Welch–Satterthwaite degrees of freedom (approximate).
    """
    nx, ny = len(x), len(y)
    mx, my = np.mean(x), np.mean(y)
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    t = (mx - my) / np.sqrt(vx / nx + vy / ny)
    df = (vx / nx + vy / ny) ** 2 / ((vx ** 2) / ((nx ** 2) * (nx - 1)) + (vy ** 2) / ((ny ** 2) * (ny - 1)))
    return t, df


def chisq_test(table):
    """
    Pearson chi-squared test on a contingency table (no p-value, no SciPy).

    Parameters
    ----------
    table : pandas.DataFrame
        Contingency table.

    Returns
    -------
    chi2 : float
        Chi-squared statistic.
    df : int
        Degrees of freedom.
    """
    obs = table.to_numpy(dtype=float)
    rsum = obs.sum(axis=1, keepdims=True)
    csum = obs.sum(axis=0, keepdims=True)
    tot = obs.sum()
    expected = rsum @ (csum / tot)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = np.nansum((obs - expected) ** 2 / expected)
    df = (obs.shape[0] - 1) * (obs.shape[1] - 1)
    return chi2, df


def pooled_sd(x, y):
    """
    Pooled standard deviation for SMD (Cohen's d).

    Parameters
    ----------
    x, y : array-like

    Returns
    -------
    float
        Pooled standard deviation.
    """
    nx, ny = len(x), len(y)
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    return np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))


def cliffs_delta(x, y):
    """
    Cliff's delta (non-parametric effect size): P(X>Y) - P(X<Y).

    Parameters
    ----------
    x, y : array-like

    Returns
    -------
    float
        Cliff's delta.
    """
    x = np.asarray(x)
    y = np.asarray(y)

    # Sort for faster counting
    x_sorted = np.sort(x)
    y_sorted = np.sort(y)

    # O(n log n) counting via binary search
    import bisect
    greater = sum(len(y_sorted) - bisect.bisect_right(y_sorted, xi) for xi in x_sorted)
    less = sum(bisect.bisect_left(y_sorted, xi) for xi in x_sorted)
    delta = (greater - less) / (len(x_sorted) * len(y_sorted))
    return delta


# ------------------------ Main ------------------------
def main():
    """
    Load and clean ABIDE phenotypic data, run basic QC summaries, and
    export a compact set of artifacts used by downstream steps.

    Outputs
    -------
    artifacts/<stage>/:
        - pheno.parquet
        - subjects.csv
        - missing_subjects.txt (if any)
        - qc_counts_by_dx.csv
        - qc_sex_by_dx.csv (if SEX available)
        - qc_site_by_dx.csv (if SITE available)
        - qc_age_by_dx_summary.csv (if AGE available)
        - qc_agebin_by_dx.csv (if AGE available)
        - QC_subgroup_variability.md (summary report)
    """
    # ------------------------ Load raw phenotype ------------------------
    df = pd.read_csv(PHENO_FILE)
    df["FILE_ID"] = df["FILE_ID"].astype(str).str.strip()

    n_before = len(df)
    df = df[df["FILE_ID"].str.lower() != "no filename"].copy()
    if len(df) < n_before:
        print(f"⚠️ {n_before - len(df)} rows removed (FILE_ID == 'no filename').")

    # ------------------------ Check presence of .1D ------------------------
    rois_dir = Path(DATA_DIR)
    available_ids = {Path(p).stem.replace("_rois_cc400", "").strip() for p in rois_dir.glob("*.1D")}
    mask = df["FILE_ID"].isin(available_ids)
    missing = df.loc[~mask, "FILE_ID"].tolist()
    df = df[mask].copy()
    if missing:
        (ART / "missing_subjects.txt").write_text("\n".join(missing))
        print(f"⚠️ {len(missing)} subjects missing .1D listed → {ART / 'missing_subjects.txt'}")

    # ------------------------ Indexing & base saves ------------------------
    df = df.set_index("FILE_ID")
    pheno_path = ART / "pheno.parquet"
    df.to_parquet(pheno_path)
    subject_path = ART / "subjects.csv"
    df[["SITE_ID", "DX_GROUP"]].reset_index().to_csv(subject_path, index=False)
    print(f"✓ cleaned phenotype → {pheno_path}")
    print(f"→ subjects.csv written → {subject_path}")

    # ------------------------ QC Preparations ------------------------
    stats_lib = try_import_scipy()
    if "DX_GROUP" not in df.columns:
        raise RuntimeError("Missing column DX_GROUP (1=ASD, 2=TDC).")

    has_age = "AGE_AT_SCAN" in df.columns
    has_sex = "SEX" in df.columns
    has_site = "SITE_ID" in df.columns

    # ------------------------ Counts by diagnosis ------------------------
    counts_dx = df["DX_GROUP"].value_counts().rename({1: "ASD", 2: "TDC"})
    counts_dx.to_csv(ART / "qc_counts_by_dx.csv")
    print("\n--- Counts by diagnosis ---")
    print(counts_dx)

    # ------------------------ SEX × DX_GROUP ------------------------
    if has_sex:
        sex_map = {1: "M", 2: "F", "M": "M", "F": "F"}
        sex = df["SEX"].map(sex_map).fillna(df["SEX"])
        tab_sex_dx = pd.crosstab(df["DX_GROUP"], sex).rename(index={1: "ASD", 2: "TDC"})
        tab_sex_dx.to_csv(ART / "qc_sex_by_dx.csv")
        print("\n--- Sex × Diagnosis ---")
        print(tab_sex_dx)

        # Proportions & alert
        props = tab_sex_dx.div(tab_sex_dx.sum(axis=1), axis=0)
        delta_m = abs(props.loc["ASD"].get("M", 0) - props.loc["TDC"].get("M", 0))
        sex_alert = delta_m >= DELTA_PROP_SEX
        print(
            f"Δ male proportion (ASD vs TDC) = {delta_m:.3f}  → "
            f"{'ALERT (imbalance)' if sex_alert else 'OK'}  (threshold {DELTA_PROP_SEX})"
        )

        # Chi-squared
        if stats_lib:
            chi2, p, dof, _ = stats_lib.chi2_contingency(tab_sex_dx.values)
            print(f"Chi²(SEX×DX) = {chi2:.2f}, dof={dof}, p={p:.3g}")
        else:
            chi2, dof = chisq_test(tab_sex_dx)
            print(f"Chi²(SEX×DX) ≈ {chi2:.2f}, dof={dof} (p not computed without SciPy)")

    # ------------------------ SITE × DX_GROUP ------------------------
    if has_site:
        tab_site_dx = pd.crosstab(df["DX_GROUP"], df["SITE_ID"]).rename(index={1: "ASD", 2: "TDC"})
        tab_site_dx.to_csv(ART / "qc_site_by_dx.csv")
        print("\n--- Site × Diagnosis (top 10 columns by total) ---")
        top_sites = tab_site_dx.sum(axis=0).sort_values(ascending=False).head(10).index
        print(tab_site_dx[top_sites])

        # Small cells
        small_cells = (tab_site_dx < N_MIN_CELL)
        n_small = int(small_cells.to_numpy().sum())
        if n_small:
            print(f"⚠️ {n_small} cells < {N_MIN_CELL} in SITE×DX (imbalance / low counts).")

        if stats_lib and tab_site_dx.shape[1] > 1:
            chi2, p, dof, _ = stats_lib.chi2_contingency(tab_site_dx.values)
            print(f"Chi²(SITE×DX) = {chi2:.2f}, dof={dof}, p={p:.3g}")
        elif tab_site_dx.shape[1] > 1:
            chi2, dof = chisq_test(tab_site_dx)
            print(f"Chi²(SITE×DX) ≈ {chi2:.2f}, dof={dof} (p not computed without SciPy)")

    # ------------------------ AGE vs DX_GROUP ------------------------
    age_summary = None
    if has_age:
        print("\n--- Age by diagnosis (summaries) ---")
        g = df.groupby("DX_GROUP")["AGE_AT_SCAN"]
        age_summary = g.agg(["count", "mean", "std", "median", "min", "max"]).rename(index={1: "ASD", 2: "TDC"})
        age_summary.to_csv(ART / "qc_age_by_dx_summary.csv")
        print(age_summary)

        # Data arrays
        x = df.loc[df["DX_GROUP"] == 1, "AGE_AT_SCAN"].dropna().values  # ASD
        y = df.loc[df["DX_GROUP"] == 2, "AGE_AT_SCAN"].dropna().values  # TDC

        if len(x) > 2 and len(y) > 2:
            # 1) Welch t-test
            if stats_lib:
                t, p_t = stats_lib.ttest_ind(x, y, equal_var=False, nan_policy="omit")
            else:
                t, _dfw = welch_ttest(x, y)
                p_t = np.nan
            # 2) Mann–Whitney U
            if stats_lib:
                u_stat, p_u = stats_lib.mannwhitneyu(x, y, alternative="two-sided")
            else:
                u_stat, p_u = np.nan, np.nan
            # 3) KS test
            if stats_lib:
                ks_stat, p_ks = stats_lib.ks_2samp(x, y, alternative="two-sided", mode="auto")
            else:
                ks_stat, p_ks = np.nan, np.nan
            # 4) Effect sizes: SMD & Cliff's delta
            sd_p = pooled_sd(x, y)
            smd = (np.mean(x) - np.mean(y)) / sd_p if sd_p > 0 else np.nan
            delta = cliffs_delta(x, y)

            print("\n--- Age tests (ASD vs TDC) ---")
            print(f"Welch t = {t:.2f}, p = {p_t:.3g}")
            print(f"Mann–Whitney U = {u_stat if not np.isnan(u_stat) else 'NA'} , p = {p_u if not np.isnan(p_u) else 'NA'}")
            print(f"KS = {ks_stat if not np.isnan(ks_stat) else 'NA'} , p = {p_ks if not np.isnan(p_ks) else 'NA'}")
            print(f"Effect sizes: SMD ≈ {smd:.2f}, Cliff’s δ ≈ {delta:.2f}")

        # 5) Age bins × DX (χ²)
        df["_AGE_BIN_"] = pd.cut(df["AGE_AT_SCAN"], bins=AGE_BINS, labels=AGE_LABELS, right=True)
        tab_agebin_dx = pd.crosstab(df["DX_GROUP"], df["_AGE_BIN_"]).rename(index={1: "ASD", 2: "TDC"})
        tab_agebin_dx.to_csv(ART / "qc_agebin_by_dx.csv")
        print("\n--- Age bins × Diagnosis (χ²) ---")
        print(tab_agebin_dx)
        if stats_lib and tab_agebin_dx.shape[1] > 1:
            chi2, p, dof, _ = stats_lib.chi2_contingency(tab_agebin_dx.fillna(0).values)
            print(f"Chi²(AgeBin×DX) = {chi2:.2f}, dof={dof}, p={p:.3g}")
        elif tab_agebin_dx.shape[1] > 1:
            chi2, dof = chisq_test(tab_agebin_dx.fillna(0))
            print(f"Chi²(AgeBin×DX) ≈ {chi2:.2f}, dof={dof} (p not computed without SciPy)")

    # ------------------------ Markdown report ------------------------
    report_md = ART / "QC_subgroup_variability.md"
    lines = []
    lines.append("# QC – Subgroup variability\n")
    lines.append(f"- Total subjects (with .1D): **{len(df)}**\n")
    lines.append(f"- Counts by diagnosis:\n\n{counts_dx.to_markdown()}\n")

    if has_sex:
        lines.append("\n## Sex × Diagnosis\n")
        lines.append(tab_sex_dx.to_markdown() + "\n")
        lines.append(
            f"- Δ male proportion (ASD vs TDC) = **{delta_m:.3f}** "
            f"→ {'ALERT' if delta_m >= DELTA_PROP_SEX else 'OK'} (threshold {DELTA_PROP_SEX})\n"
        )

    if has_site:
        lines.append("\n## Site × Diagnosis (top 10)\n")
        lines.append(tab_site_dx[top_sites].to_markdown() + "\n")
        if n_small:
            lines.append(f"- ⚠️ **{n_small}** cells < **{N_MIN_CELL}** in SITE×DX (imbalance / low counts).\n")

    if has_age:
        lines.append("\n## Age by diagnosis – summaries\n")
        lines.append(age_summary.to_markdown() + "\n")
        lines.append("\n## Age – tests\n")
        if len(x) > 2 and len(y) > 2:
            lines.append(f"- Welch t = **{t:.2f}**, p = **{p_t:.3g}**\n")
            lines.append(f"- Mann–Whitney U = **{u_stat if not np.isnan(u_stat) else 'NA'}**, p = **{p_u if not np.isnan(p_u) else 'NA'}**\n")
            lines.append(f"- KS = **{ks_stat if not np.isnan(ks_stat) else 'NA'}**, p = **{p_ks if not np.isnan(p_ks) else 'NA'}**\n")
            lines.append(f"- Effect sizes: SMD ≈ **{smd:.2f}**, Cliff’s δ ≈ **{delta:.2f}**\n")
        lines.append("\n## Age bins × Diagnosis (χ²)\n")
        lines.append(tab_agebin_dx.to_markdown() + "\n")

    report_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"\n📝 QC report written → {report_md}")


# ------------------------ Entrypoint ------------------------
if __name__ == "__main__":
    main()
