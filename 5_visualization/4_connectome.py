# 5_visualization_generic/plot_connectome_network.py
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import numpy as np
from pathlib import Path
from nilearn import plotting
import matplotlib.cm as cm

# ------------------------ Paths / Folders ------------------------
ART_VIZ = Path("artifacts/5_visualization")
CONNECTOME_DIR = ART_VIZ / "connectome"
CONNECTOME_DIR.mkdir(parents=True, exist_ok=True)

# -- Inputs --
NETWORKS_CSV = Path("artifacts/4_statistics/network_tests.csv")
ROI_TO_YEO = Path("artifacts/3_annotation_networks/roi_to_yeo.csv")
CRADD_LABELS = Path("atlases/atlas_craddock400/CC400_ROI_labels.csv")

# ------------------------ 1) Load network_tests ------------------------
df_net = pd.read_csv(NETWORKS_CSV)

# Normalize a few column names (spaces, dashes)
df_net = df_net.rename(columns={
    "t value": "t_value",
    "p-value": "p_value",
    "p_FDR": "p_FDR",
    "Connection": "Connection",
    "type": "type",
})

# Basic cleaning
df_net["type"] = df_net["type"].str.strip().str.lower()
df_net["Connection"] = df_net["Connection"].str.strip()

# Abbreviation -> full Yeo-7 names used for the figure
abbr2full = {
    "VIS": "Visual",
    "SMN": "Somatomotor",
    "DAN": "Dorsal Attention",
    "SAL": "Salience / Ventral Attention",
    "LIM": "Limbic",
    "FPN": "Frontoparietal Control",
    "DMN": "Default Mode",
}
nets = [
    "Visual",
    "Somatomotor",
    "Dorsal Attention",
    "Salience / Ventral Attention",
    "Limbic",
    "Frontoparietal Control",
    "Default Mode",
]

# ------------------------ 2) Build M (network×network) from pairs ------------------------

# >>> Optional filtering <<<
# Set this flag to True if you want only significant connections
ONLY_SIGNIF = True

if ONLY_SIGNIF:
    df_used = df_net[df_net["p_FDR"] < 0.05]
    print(f"[INFO] Retained connections: {len(df_used)} / {len(df_net)} (p_FDR < 0.05)")
else:
    df_used = df_net.copy()

# Initialize matrix with 0.0
M_df = pd.DataFrame(0.0, index=nets, columns=nets, dtype=float)

for _, row in df_used.iterrows():
    tval = float(row["t_value"])
    if row["type"] == "inter":
        try:
            a, b = [x.strip() for x in row["Connection"].split("-")]
        except Exception:
            continue
        if a in abbr2full and b in abbr2full:
            A, B = abbr2full[a], abbr2full[b]
            M_df.loc[A, B] = tval
            M_df.loc[B, A] = tval
    elif row["type"] == "intra":
        a = row["Connection"].strip()
        if a in abbr2full:
            A = abbr2full[a]
            M_df.loc[A, A] = tval

M = M_df.values.astype(float)

# ------------------------ 3) Node color = intra-network t (diagonal) ------------------------
node_color_vals = np.diag(M).copy()

# ------------------------ 4) Network coords = median of ROI centers ------------------------
map_roi = pd.read_csv(ROI_TO_YEO).set_index("ROI_number")
labs = pd.read_csv(CRADD_LABELS)

def parse_com(s: str):
    """Extract (x,y,z) from a string like '(28.6;55.2;18.2)' or '(28.6,55.2,18.2)'; NaN on failure."""
    import re
    if not isinstance(s, str):
        return np.nan, np.nan, np.nan
    nums = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', s.replace(';', ','))
    if len(nums) < 3:
        return np.nan, np.nan, np.nan
    x, y, z = map(float, nums[:3])
    return x, y, z

xyz = labs["center_of_mass"].apply(parse_com).apply(pd.Series)
xyz.columns = ["x", "y", "z"]
labs = pd.concat([labs[["ROI_number"]], xyz], axis=1).rename(columns={"ROI_number": "roi_number"})
labs = labs.merge(map_roi[["Yeo7_id"]], left_on="roi_number", right_index=True, how="left")
labs["Yeo7_id"] = labs["Yeo7_id"].fillna(0).astype(int)

yeo7_names = {
    1: "Visual",
    2: "Somatomotor",
    3: "Dorsal Attention",
    4: "Salience / Ventral Attention",
    5: "Limbic",
    6: "Frontoparietal Control",
    7: "Default Mode",
}

coords = []
for k in range(1, 8):
    xyz_k = labs.query("Yeo7_id == @k")[['x', 'y', 'z']].dropna()
    coords.append(tuple(np.median(xyz_k.values, axis=0)) if len(xyz_k) else (np.nan, np.nan, np.nan))

# Order as 'nets'
order_names = nets  # same order as M_df
# Build coords in the same order as M_df
name2coord = {yeo7_names[i + 1]: coords[i] for i in range(7)}
coords = [name2coord[nm] for nm in order_names]

# ------------------------ 5) Adjacency without self-loops ------------------------
adj = M.copy()
np.fill_diagonal(adj, 0.0)

coords_arr = np.array(coords, dtype=float)
n_valid = np.isfinite(coords_arr).all(axis=1).sum()
print(f"[connectome] networks with valid coordinates: {n_valid}/7")
if n_valid < 3:
    raise RuntimeError(
        "Too few networks have valid coordinates. Check CRADD_LABELS.center_of_mass and parsing."
    )

# ------------------------ 6) Plot & save ------------------------
out_png = CONNECTOME_DIR / "fig_connectome_yeo7.png"
try:
    # Palette (7 distinct colors)
    cmap_nodes = cm.get_cmap("Set2", 7)
    node_colors = [cmap_nodes(i) for i in range(7)]  # same order as 'nets'

    plot = plotting.plot_connectome(
        adj, coords,
        node_color=node_colors,
        node_size=60,
        edge_cmap="coolwarm",
        edge_vmin=-np.nanmax(np.abs(adj)) if np.isfinite(adj).any() else -1,
        edge_vmax= np.nanmax(np.abs(adj)) if np.isfinite(adj).any() else  1,
        colorbar=True,
        title="Network connectome (Yeo-7) — edges = t (ASD − TDC)",
    )

    # Legend
    handles = [mpatches.Patch(color=node_colors[i], label=nets[i]) for i in range(7)]
    fig = mpl.pyplot.gcf()
    fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(-0.04, -0.03),
        ncol=3,
        frameon=False,
    )

    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"✓ Connectome saved → {out_png}")
except Exception as e:
    print(f"[ERROR] connectome save: {e}")
    raise
