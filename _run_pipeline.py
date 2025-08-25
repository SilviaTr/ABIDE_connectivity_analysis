# _run_pipeline.py
import os
import subprocess
import sys
from pathlib import Path

# repo root with config.py
from config import ROOT  # ROOT should be a Path pointing to the repo root

def run_script(script_rel_path: str, *args: str) -> None:
    """
    Run a pipeline step as a subprocess while ensuring the repo root is on PYTHONPATH.
    Example: run_script("1_data_cleaning/1_load_and_clean_pheno.py")
    """
    script_path = (ROOT / script_rel_path).resolve()
    if not script_path.exists():
        raise FileNotFoundError(f"Script not found: {script_path}")

    cmd = [sys.executable, str(script_path), *map(str, args)]

    # Ensure imports like `from config import ...` work in every step
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{ROOT}{os.pathsep}{env.get('PYTHONPATH', '')}"

    print(f"\n=== RUN: {' '.join(cmd)} ===")
    # Set cwd to ROOT for good measure (some scripts may use relative paths)
    subprocess.run(cmd, check=True, cwd=str(ROOT), env=env)

def main():
    # 0) Download preprocessed data (defaults are inside the script)
    run_script("0_download_data/download_abide_preprocessed_data.py")

    # 1) Phenotypic cleaning
    run_script("1_data_cleaning/1_load_and_clean_pheno.py")

    # 2) Connectivity extraction
    run_script("2_connectivity_extraction/1_compute_connectivity.py")
    run_script("2_connectivity_extraction/2_regress_covariates.py")

    # 3) ROI / network annotations
    run_script("3_annotation_networks/1_map_roi_to_HO.py")
    run_script("3_annotation_networks/2_map_roi_to_yeo.py")

    # 4) Statistics
    run_script(
        "4_statistics/1_build_network_scores.py",
        "--sparsity", "0.2",
        "--n-pca", "3",
        "--pca-report",
        "--anchor-mode", "signed"  # explicit for reproducibility
    )
    run_script("4_statistics/2_run_ttest.py")
    # Optional: edge-level stats
    # run_script("4_statistics/compute_edge_ttests_fdr.py")

    # 5) Visualization
    run_script("5_visualization/1_summary_table.py")
    run_script("5_visualization/2_build_heatmap.py")
    run_script("5_visualization/3_1_compute_tscore_roi.py")
    run_script("5_visualization/3_2_plot_surface_tscore.py")
    run_script("5_visualization/4_connectome.py")

    print("\n✓ Full ABIDE pipeline executed successfully.")

if __name__ == "__main__":
    main()
