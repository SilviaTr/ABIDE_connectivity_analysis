# 5_visualization/1_summary_table.py

import pandas as pd
from pathlib import Path

# ------------------------ Folders ------------------------
ART_VIZ = Path("artifacts/5_visualization")
TABLE_DIR = ART_VIZ / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

IN_CSV = Path("artifacts/4_statistics/network_tests.csv")  # intra/inter blocks

# ------------------------ Threshold ------------------------
ALPHA = 0.05  # threshold on p_FDR


# ------------------------ Helpers ------------------------
def latex_escape(s: str) -> str:
    """Escape LaTeX special characters in a string."""
    repl = {
        "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
        "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}", "\\": r"\textbackslash{}",
    }
    return "".join(repl.get(ch, ch) for ch in str(s))


def fmt_t(v):
    """Format t value with 2 decimals if numeric."""
    try:
        return f"{float(v):.2f}"
    except Exception:
        return str(v)


def tabular_from_rows(rows, caption, label):
    """Build a LaTeX table environment from rows."""
    lines = []
    lines.append(r"\begin{table}[h!]")
    lines.append(r"\centering")
    lines.append(rf"\caption{{{latex_escape(caption)}}}")
    lines.append(rf"\label{{{latex_escape(label)}}}")
    lines.append(r"\begin{tabular}{lcccc}")
    lines.append(r"\toprule")
    lines.append(r"Connection & HC & ASD & t value & p-value \\")
    lines.append(r"\midrule")
    for r in rows:
        lines.append(
            f"{latex_escape(r['Connection'])} & {latex_escape(r['HC'])} & "
            f"{latex_escape(r['ASD'])} & {latex_escape(r['t value'])} & "
            f"{latex_escape(r['p-value'])} \\\\"
        )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def build_table(df, block_type: str, caption: str, label: str):
    """Filter by block type, keep FDR-significant rows, and render a LaTeX table."""
    sub = df[df["type"] == block_type].copy()
    sub["p_FDR"] = pd.to_numeric(sub["p_FDR"], errors="coerce")
    sub = sub[sub["p_FDR"] < ALPHA].sort_values("Connection")

    rows = []
    for _, r in sub.iterrows():
        rows.append({
            "Connection": r["Connection"],
            "HC": r["HC"],
            "ASD": r["ASD"],
            "t value": fmt_t(r["t value"]),
            "p-value": r.get("p-value", ""),
        })
    return tabular_from_rows(rows, caption, label)


def assemble_document(intra_block: str, inter_block: str) -> str:
    """Assemble a minimal English LaTeX document with two tables."""
    preamble = r"""\documentclass[a4paper,11pt]{article}
\usepackage{booktabs}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[english]{babel}
\begin{document}
\section*{Network-level connectivity tables}
"""
    enddoc = r"\end{document}"
    return "\n".join([preamble, intra_block, "", inter_block, enddoc])


# ------------------------ Main ------------------------
def main():
    assert IN_CSV.exists(), f"Not found: {IN_CSV}"
    df = pd.read_csv(IN_CSV)

    # numeric conversions (if present)
    for col in ("t value", "p_FDR"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    intra_block = build_table(df, "intra", "Intra-network connectivity (Yeo-7)", "tab:intra_blocks")
    inter_block = build_table(df, "inter", "Inter-network connectivity (Yeo-7)", "tab:inter_blocks")

    out_path = TABLE_DIR / "connectivity_tables_blocks.tex"
    out_path.write_text(assemble_document(intra_block, inter_block), encoding="utf-8")
    print(f"✓ LaTeX written → {out_path}")


if __name__ == "__main__":
    main()
