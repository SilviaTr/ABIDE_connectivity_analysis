# 5_visualization_generic/2_build_heatmap.py

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ------------------------ Folders ------------------------
ART_VIZ = Path("artifacts/5_visualization")
HEATMAP_DIR = ART_VIZ / "heatmap"
HEATMAP_DIR.mkdir(parents=True, exist_ok=True)

IN_CSV = Path("artifacts/4_statistics/network_tests.csv")
ALPHA = 0.05

# ------------------------ Network names ------------------------
ABBR2FULL = {
    "VIS": "Visual",
    "SMN": "Somatomotor",
    "DAN": "Dorsal Attention",
    "SAL": "Salience / Ventral Attention",
    "LIM": "Limbic",
    "FPN": "Frontoparietal Control",
    "DMN": "Default Mode",
}
YEO7_ORDER = [
    "Visual", "Somatomotor", "Dorsal Attention",
    "Salience / Ventral Attention", "Limbic",
    "Frontoparietal Control", "Default Mode",
]

# Placeholder labels sometimes present in Yeo-17 mapping
EXCLUDE_NETS = {"A", "B", "C"}


# ------------------------ Helpers ------------------------
def expand_name(x: str) -> str:
    """Expand abbreviations (e.g., DMN → Default Mode)."""
    x = ("" if pd.isna(x) else str(x)).strip()
    return ABBR2FULL.get(x, x)


def split_pair(s: str):
    """Split a connection 'A-B' into (A,B). For intra, return (A,None)."""
    s = ("" if pd.isna(s) else str(s)).strip()
    for sep in ["—", " – ", "-", " –", "– ", " - ", " -", "- "]:
        if sep in s:
            a, b = s.split(sep, 1)
            return a.strip(), b.strip()
    return s, None


def detect_parcellation(nets_present):
    """Detect if the set corresponds to Yeo-7 or Yeo-17."""
    nets_present = {expand_name(n) for n in nets_present if str(n).strip() != ""}
    if nets_present and nets_present.issubset(set(YEO7_ORDER)):
        return "yeo7", YEO7_ORDER[:]
    nets17 = sorted(n for n in nets_present if n not in EXCLUDE_NETS)
    return "yeo17", nets17


def ensure_square(order):
    """Initialize a symmetric DataFrame filled with NaN."""
    return pd.DataFrame(np.nan, index=order, columns=order, dtype=float)


# ------------------------ Main ------------------------
def main():
    assert IN_CSV.exists(), f"Not found: {IN_CSV}"
    df = pd.read_csv(IN_CSV)

    # Ensure numeric
    for col in ["t value", "p_FDR"]:
        if col not in df.columns:
            raise ValueError(f"Missing column: {col}")
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # If type missing, infer from connection string
    if "type" not in df.columns:
        df["type"] = df["Connection"].astype(str).apply(
            lambda c: "inter" if split_pair(c)[1] else "intra"
        )

    # Expand network names
    A, B, nets_seen = [], [], set()
    for conn in df["Connection"].astype(str):
        a_raw, b_raw = split_pair(conn)
        a = expand_name(a_raw)
        b = expand_name(b_raw) if b_raw is not None else None
        A.append(a)
        B.append(b)
        if a:
            nets_seen.add(a)
        if b:
            nets_seen.add(b)
    df["A"], df["B"] = A, B

    # Detect parcellation & order
    yeo_version, ORDER = detect_parcellation(nets_seen)
    print(f"→ Detected parcellation: {yeo_version.upper()} | #networks = {len(ORDER)}")

    # Exclude placeholders if Yeo-17
    if yeo_version == "yeo17":
        before = len(df)
        df = df[
            ~df["A"].isin(EXCLUDE_NETS) & (~df["B"].isin(EXCLUDE_NETS) | df["B"].isna())
        ].copy()
        removed = before - len(df)
        if removed > 0:
            print(f"  ↳ Removed {removed} rows with placeholders {sorted(EXCLUDE_NETS)}")

    # Square matrix
    M = ensure_square(ORDER)

    # Fill intra-network (significant only)
    intra = df[(df["type"] == "intra") & (df["p_FDR"] < ALPHA)].copy()
    for _, r in intra.iterrows():
        a = r["A"]
        if a in M.index:
            M.loc[a, a] = float(r["t value"])

    # Fill inter-network (significant only)
    inter = df[(df["type"] == "inter") & (df["p_FDR"] < ALPHA)].copy()
    for _, r in inter.iterrows():
        a, b = r["A"], r["B"]
        if a in M.index and b in M.columns:
            val = float(r["t value"])
            M.loc[a, b] = val
            M.loc[b, a] = val

    M_ord = M.loc[ORDER, ORDER]

    # Plot heatmap
    plt.figure(figsize=(max(6, 0.7 * len(ORDER)), max(5.5, 0.7 * len(ORDER))))
    mask = ~np.isfinite(M_ord.values)
    sns.heatmap(
        M_ord,
        cmap="coolwarm",
        center=0,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "t (ASD − HC)"},
        mask=mask,
    )
    plt.title(f"Significant blocks (p_FDR < {ALPHA}) — t (ASD−HC) — {yeo_version.upper()}")
    plt.tight_layout()

    out_png = HEATMAP_DIR / f"heatmap_tvalues_blocks_{yeo_version}_sig.png"
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"✓ Heatmap saved → {out_png}")


if __name__ == "__main__":
    main()
