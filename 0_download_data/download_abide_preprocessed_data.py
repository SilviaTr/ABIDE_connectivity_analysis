# 0_download_data/download_abide_preprocessed_data.py
"""
Download ABIDE preprocessed derivatives (e.g., CPAC rois_cc400) with motion-QC filtering.

Defaults:
- pipeline:  cpac
- strategy:  filt_global
- derivative: rois_cc400  (ROI time series, .1D)
- diagnosis filter: both  (ASD + TDC)
- output dir: <repo_root>/0_download_data/preprocessed_dataset

Usage examples:
    python 0_download_data/download_abide_preprocessed_data.py
    python 0_download_data/download_abide_preprocessed_data.py \
        --pipeline cpac --strategy filt_global --derivative rois_cc400 --diagnosis both
    python 0_download_data/download_abide_preprocessed_data.py --out-dir /data/ABIDE/preproc
"""

from __future__ import annotations
import argparse
import os
from pathlib import Path
import urllib.request as request

try:
    # Prefer using your config if available (ROOT points at the repo root)
    from config import ROOT as REPO_ROOT
except Exception:
    # Fallback: repo root assumed to be the parent of this file's parent (…/0_download_data/this_file.py)
    REPO_ROOT = Path(__file__).resolve().parents[1]


S3_PREFIX = "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative"
S3_PHENO  = f"{S3_PREFIX}/Phenotypic_V1_0b_preprocessed1.csv"
DEFAULT_OUT = REPO_ROOT / "0_download_data" / "preprocessed_dataset"
LOCAL_PHENO = REPO_ROOT / "0_download_data" / "Phenotypic_V1_0b_preprocessed1.csv"
MEAN_FD_THRESH = 0.2  # motion QC


def ensure_out_dir(out_dir: Path) -> Path:
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def safe_remove_prefix(s: str, prefix: str) -> str:
    return s[len(prefix):] if s.startswith(prefix) else s


def download_text(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with request.urlopen(url) as r, open(dest, "wb") as f:
        f.write(r.read())


def collect_and_download(
    derivative: str,
    pipeline: str,
    strategy: str,
    out_dir: Path,
    diagnosis: str,
) -> None:
    """
    Collect S3 paths for ABIDE preprocessed derivatives and download them locally.

    Parameters
    ----------
    derivative : str
        e.g., 'rois_cc400', 'reho', …
        If it contains 'roi', files are expected to end with '.1D', else '.nii.gz'.
    pipeline : str
        e.g., 'cpac', 'ccs', 'dparsf', 'niak'
    strategy : str
        e.g., 'filt_global', 'filt_noglobal', 'nofilt_global', 'nofilt_noglobal'
    out_dir : Path
        Local destination directory (will be created).
    diagnosis : str
        One of {'asd', 'tdc', 'both'} — filters by DX_GROUP (1=ASD, 2=TDC).
    """
    derivative = derivative.lower()
    pipeline = pipeline.lower()
    strategy = strategy.lower()
    diagnosis = diagnosis.lower()

    # Decide extension by derivative type
    extension = ".1D" if "roi" in derivative else ".nii.gz"

    out_dir = ensure_out_dir(out_dir)

    # Get phenotypic CSV (cache a local copy under 0_download_data/)
    if not LOCAL_PHENO.exists():
        print(f"[info] Downloading phenotypic CSV → {LOCAL_PHENO}")
        download_text(S3_PHENO, LOCAL_PHENO)
    else:
        print(f"[info] Using cached phenotypic CSV → {LOCAL_PHENO}")

    # Read phenotypic CSV lines
    with open(LOCAL_PHENO, "rb") as f:
        pheno_lines = f.readlines()

    header = pheno_lines[0].decode().strip().split(",")
    try:
        site_idx = header.index("SITE_ID")
        file_idx = header.index("FILE_ID")
        age_idx = header.index("AGE_AT_SCAN")
        sex_idx = header.index("SEX")
        dx_idx = header.index("DX_GROUP")
        mean_fd_idx = header.index("func_mean_fd")
    except Exception as e:
        raise RuntimeError("Failed to parse phenotypic header") from e

    # Build S3 paths with filtering
    s3_paths = []
    for row in pheno_lines[1:]:
        cs = row.decode().strip().split(",")

        try:
            row_file_id = cs[file_idx]
            row_mean_fd = float(cs[mean_fd_idx])
            row_dx = cs[dx_idx]
        except Exception:
            continue  # skip malformed rows

        if row_file_id == "no_filename":
            continue
        if row_mean_fd >= MEAN_FD_THRESH:
            continue

        if (diagnosis == "asd" and row_dx != "1") or (diagnosis == "tdc" and row_dx != "2"):
            continue

        filename = f"{row_file_id}_{derivative}{extension}"
        s3 = "/".join([S3_PREFIX, "Outputs", pipeline, strategy, derivative, filename])
        s3_paths.append(s3)

    # Download loop
    total = len(s3_paths)
    if total == 0:
        print("[warn] No files matched the given filters. Nothing to download.")
        return

    for k, s3_url in enumerate(s3_paths, start=1):
        # Keep the S3 path sub-structure under out_dir
        # e.g., Outputs/cpac/filt_global/rois_cc400/subject_rois_cc400.1D
        rel = safe_remove_prefix(s3_url, S3_PREFIX)
        rel = rel.lstrip("/")  # ensure no leading slash
        dst = out_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)

        if dst.exists():
            print(f"[skip] {dst}")
            continue

        try:
            print(f"[{k:5d}/{total:5d}] Downloading → {dst}")
            request.urlretrieve(s3_url, dst)
        except Exception as e:
            print(f"[error] Failed: {s3_url} ({e})")

    print("Done.")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Download ABIDE preprocessed derivatives.")
    p.add_argument("--pipeline", default="cpac",
                   choices=["cpac", "ccs", "dparsf", "niak"],
                   help="Preprocessing pipeline.")
    p.add_argument("--strategy", default="filt_global",
                   choices=["filt_global", "filt_noglobal", "nofilt_global", "nofilt_noglobal"],
                   help="Denoising strategy.")
    p.add_argument("--derivative", default="rois_cc400",
                   help="Derivative name (e.g., rois_cc400, reho, vmhc).")
    p.add_argument("--diagnosis", default="both", choices=["asd", "tdc", "both"],
                   help="Filter by diagnosis: ASD, TDC or both.")
    p.add_argument("--out-dir", default=str(DEFAULT_OUT),
                   help="Destination directory. Defaults to 0_download_data/preprocessed_dataset under repo root.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Resolve output directory: if user passes a relative path (e.g., "preprocessed_dataset"),
    # place it under 0_download_data/ relative to the repo root.
    out_dir_arg = Path(args.out_dir)
    if not out_dir_arg.is_absolute():
        # If it's exactly the default string, DEFAULT_OUT already points to the correct place.
        # Otherwise, put it under 0_download_data/<given_relative>
        if out_dir_arg == Path("preprocessed_dataset"):
            out_dir = DEFAULT_OUT
        else:
            out_dir = (REPO_ROOT / "0_download_data" / out_dir_arg).resolve()
    else:
        out_dir = out_dir_arg

    collect_and_download(
        derivative=args.derivative,
        pipeline=args.pipeline,
        strategy=args.strategy,
        out_dir=out_dir,
        diagnosis=args.diagnosis,
    )
