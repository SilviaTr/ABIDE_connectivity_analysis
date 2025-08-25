# 3_annotation_networks/2_map_roi_to_yeo.py

import json
from pathlib import Path

import numpy as np
import pandas as pd
from nilearn import datasets
from nilearn.image import load_img, resample_to_img, new_img_like

from config import ATLAS_PATH, YEO7_PATH, YEO17_PATH, stage_artifacts

# ------------------------ Artifacts & Folders ------------------------
ART_CONN = stage_artifacts("2_connectivity_extraction")
CONN_DIR = ART_CONN / "connectivity"
ART_ANN = stage_artifacts("3_annotation_networks")
ART_ANN.mkdir(parents=True, exist_ok=True)

# ------------------------ Canonical Yeo network names ------------------------
YEO7_NAMES = {
    1: "Visual",
    2: "Somatomotor",
    3: "Dorsal Attention",
    4: "Salience / Ventral Attention",
    5: "Limbic",
    6: "Frontoparietal Control",
    7: "Default Mode",
}
YEO17_NAMES = {
     1: "Visual – Central",    2: "Visual – Peripheral",
     3: "Somatomotor – A",     4: "Somatomotor – B",
     5: "Dorsal Attention – A", 6: "Dorsal Attention – B",
     7: "Salience/Ventral Attention – A", 8: "Salience/Ventral Attention – B",
     9: "Limbic – A",         10: "Limbic – B",
    11: "Frontoparietal Control – A", 12: "Frontoparietal Control – B", 13: "Frontoparietal Control – C",
    14: "Default Mode – A",   15: "Default Mode – B", 16: "Default Mode – C",
    17: "Temporo-Parietal",
}


# ------------------------ Helpers ------------------------
def _get_src_img(path: Path, which: str):
    """
    Load a Yeo atlas image, using local path if available, otherwise nilearn fallback.

    Parameters
    ----------
    path : Path
        Expected local path to Yeo atlas (7 or 17).
    which : str
        "7" or "17" to pick the proper nilearn atlas if local file is missing.

    Returns
    -------
    Nifti1Image
        Source Yeo atlas image.
    """
    if path.exists():
        return load_img(str(path))
    yeo = datasets.fetch_atlas_yeo_2011()
    return load_img(yeo.thick_7 if which == "7" else yeo.thick_17)


def _resample_to_atlas(src_img, atlas_img, dst_path: Path, max_label: int):
    """
    Resample a source atlas to the reference atlas grid with nearest-neighbor,
    sanitize labels, and save to disk.

    Parameters
    ----------
    src_img : Nifti1Image
        Source image to resample.
    atlas_img : Nifti1Image
        Target/reference atlas image (voxel grid + affine).
    dst_path : Path
        Output path where the resampled image will be written.
    max_label : int
        Maximum valid label value (inclusive).

    Returns
    -------
    Nifti1Image
        Resampled image loaded back from disk (for consistency).
    """
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    res = resample_to_img(
        src_img,
        atlas_img,
        interpolation="nearest",
        force_resample=True,
        copy_header=True,
    )
    data = np.asanyarray(res.get_fdata())
    if data.ndim == 4 and data.shape[-1] == 1:
        data = data[..., 0]
    data = np.rint(data).astype(np.int16)
    data[data < 0] = 0
    data = np.clip(data, 0, max_label)
    out = new_img_like(atlas_img, data, copy_header=True)
    out.header.set_data_dtype(np.int16)
    out.to_filename(str(dst_path))
    return load_img(str(dst_path))


def _majority_label(parc_int: np.ndarray, mask: np.ndarray) -> int:
    """
    Return the most frequent (non-zero) label in `parc_int` within `mask`.

    Parameters
    ----------
    parc_int : np.ndarray
        Integer-labeled parcellation volume.
    mask : np.ndarray (bool)
        Boolean mask selecting voxels that belong to a given ROI.

    Returns
    -------
    int
        Majority label (0 if no voxels labeled).
    """
    vals = parc_int[mask]
    vals = vals[vals > 0]
    if vals.size == 0:
        return 0
    ids, cnt = np.unique(vals, return_counts=True)
    return int(ids[np.argmax(cnt)])


# ------------------------ Main ------------------------
def main():
    """
    Build ROI → Yeo7/Yeo17 network mapping by majority vote within each ROI.

    Steps
    -----
    1) Read n_rois from connectivity metadata.
    2) Load reference atlas (Craddock or provided) and Yeo 7/17 atlases, resampled to the reference grid.
    3) For each ROI number in 1..n_rois, compute majority Yeo7 and Yeo17 labels.
    4) Save a CSV mapping (roi_to_yeo.csv) with ROI_number, Yeo7_id/name, Yeo17_id/name.
    """
    # ------------------------ Read n_rois from connectivity metadata ------------------------
    shape = json.loads((CONN_DIR / "conn_shape.json").read_text())["shape"]
    n_rois = int(shape[1])

    # ------------------------ Load reference and Yeo atlases ------------------------
    atlas_img = load_img(str(ATLAS_PATH))
    crad = atlas_img.get_fdata().astype(int)

    y7_img = _resample_to_atlas(_get_src_img(Path(YEO7_PATH), "7"), atlas_img, Path(YEO7_PATH), 7)
    y17_img = _resample_to_atlas(_get_src_img(Path(YEO17_PATH), "17"), atlas_img, Path(YEO17_PATH), 17)

    assert y7_img.shape == atlas_img.shape and y17_img.shape == atlas_img.shape, "Resampled Yeo atlases must match reference shape."

    y7 = y7_img.get_fdata().astype(int)
    y17 = y17_img.get_fdata().astype(int)

    # ------------------------ Majority vote per ROI ------------------------
    rows = []
    for roi in range(1, n_rois + 1):
        mask = (crad == roi)
        y7_id = _majority_label(y7, mask)
        y17_id = _majority_label(y17, mask)
        rows.append((
            roi,
            y7_id,
            YEO7_NAMES.get(y7_id, "Unknown") if y7_id else "Unknown",
            y17_id,
            YEO17_NAMES.get(y17_id, "Unknown") if y17_id else "Unknown",
        ))

    df = pd.DataFrame(rows, columns=["ROI_number", "Yeo7_id", "Yeo7_name", "Yeo17_id", "Yeo17_name"])

    # ------------------------ Save mapping ------------------------
    out = ART_ANN / "roi_to_yeo.csv"
    df.to_csv(out, index=False)
    print(f"✓ ROI → Yeo7/Yeo17 mapping written → {out}")


# ------------------------ Entrypoint ------------------------
if __name__ == "__main__":
    main()
