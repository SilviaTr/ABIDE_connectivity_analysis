# config.py
from pathlib import Path

# Project root = folder that contains this config.py
ROOT = Path(__file__).resolve().parent

# ---- INPUTS (from 0_download_data) ----
DATA_DIR   = ROOT / "0_download_data" / "preprocessed_dataset" / "Outputs" / "cpac" / "filt_global" / "rois_cc400"
PHENO_FILE = ROOT / "0_download_data" / "Phenotypic_V1_0b_preprocessed1.csv"
CONN_CSV_FILE = ROOT / "artifacts" / "4_annotation_networks"/"df_conn_annotated.csv"
CONN_CSV_FILE_BOTH = ROOT / "artifacts" / "4_annotation_networks"/ "df_conn_annotated.csv"
TIMESERIES_DIR = ROOT / "0_download_data" / "preprocessed_dataset" / "Outputs" / "cpac" / "filt_global" / "rois_cc400"
TIMESERIES_DIR.mkdir(parents=True, exist_ok=True)

# ---- ATLAS & LABELS (adjust if you placed them elsewhere) ----
ATLAS_PATH = ROOT / "atlases" / "atlas_craddock400" / "cc400_roi_atlas.nii.gz"
LABELS_CSV = ROOT / "atlases" / "atlas_craddock400" / "CC400_ROI_labels.csv"

YEO7_PATH  = ROOT / "atlases" / "yeo" / "yeo7_on_craddock.nii.gz"
YEO17_PATH = ROOT / "atlases" / "yeo" / "yeo17_on_craddock.nii.gz"


# ---- OUTPUTS ----
ARTIFACTS = ROOT / "artifacts"; ARTIFACTS.mkdir(exist_ok=True, parents=True)

# ---- STATS ----
ALPHA = 0.05
FDR_METHOD = "fdr_bh"
USE_EQUAL_VAR = False
SEED = 42
COHEN_D_MIN = 0.3


# ---- VISU ----
TOP_N_EDGES = 10
SURF_FSAVERAGE = "fsaverage"

# ---- SEUIL ----
MAX_BAD_ROIS_RATIO = 0.15   # exclure si > 15% des ROIs sont "bad"

def stage_artifacts(stage_name: str) -> Path:
    """Retourne le dossier artifacts d'une étape et le crée au besoin."""
    p = ARTIFACTS / stage_name
    p.mkdir(parents=True, exist_ok=True)
    return p
