#!/usr/bin/env python3
"""Generate article-style figures for the PRJNA400142 poultry 16S pipeline."""

from __future__ import annotations

import io
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal

try:
    from matplotlib_venn import venn2
except Exception:  # pragma: no cover
    venn2 = None

PROJECT_DIR = Path(__file__).resolve().parents[1]
EXPORTED_DIR = PROJECT_DIR / "results" / "exported"
FIGURES_DIR = PROJECT_DIR / "results" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

TAXA_L3 = EXPORTED_DIR / "taxa_level3.tsv"
TAXA_L6 = EXPORTED_DIR / "taxa_level6.tsv"
METADATA = EXPORTED_DIR / "sample-metadata.used.tsv"
SHANNON = EXPORTED_DIR / "shannon" / "alpha-diversity.tsv"


def read_biom_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    kept = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("# Constructed"):
                continue
            if line.startswith("#OTU ID"):
                line = line.replace("#OTU ID", "Taxon", 1)
            kept.append(line)
    df = pd.read_csv(io.StringIO("".join(kept)), sep="\t")
    if "Taxon" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Taxon"})
    for col in df.columns[1:]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    return df


def read_metadata(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    lines = []
    with path.open() as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                lines.append(line)
    md = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str).fillna("")
    if "sample-id" not in md.columns:
        md = md.rename(columns={md.columns[0]: "sample-id"})
    if "body_site" not in md.columns:
        raise ValueError("Metadata must contain a body_site column with Lung or Trachea.")
    md["body_site"] = md["body_site"].map(normalize_site)
    md = md[md["body_site"].isin(["Lung", "Trachea"])]
    if md.empty:
        raise ValueError("No Lung/Trachea samples found in metadata.")
    return md


def normalize_site(value: str) -> str:
    text = str(value).strip().lower()
    if text in {"lung", "lungs", "pulmao", "pulmão"}:
        return "Lung"
    if text in {"trachea", "tracheal", "traqueia"}:
        return "Trachea"
    return str(value).strip()


def clean_taxon_label(taxon: str) -> str:
    bad = {"", "__", "unassigned", "unclassified", "uncultured", "unknown", "nan", "none"}
    parts = [p.strip() for p in str(taxon).split(";")]
    for part in reversed(parts):
        label = re.sub(r"^[a-zA-Z]__", "", part).strip()
        if label.lower() not in bad:
            return label
    return "Unclassified"


def group_by_site(table: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    sample_cols = [c for c in table.columns if c != "Taxon"]
    md = metadata.set_index("sample-id")
    used_cols = [c for c in sample_cols if c in md.index]
    if not used_cols:
        raise ValueError("No table sample columns match metadata sample-id values.")
    x = table[["Taxon"] + used_cols].copy()
    x["Taxon"] = x["Taxon"].map(clean_taxon_label)
    x = x.groupby("Taxon", as_index=True)[used_cols].sum()
    out = pd.DataFrame(index=x.index)
    for site in ["Lung", "Trachea"]:
        cols = [c for c in used_cols if md.loc[c, "body_site"] == site]
        out[site] = x[cols].sum(axis=1) if cols else 0
    out = out.groupby(out.index).sum()
    out = out.loc[out.sum(axis=1).sort_values(ascending=False).index]
    return out


def relative_top_table(counts: pd.DataFrame, top_n: int) -> pd.DataFrame:
    rel = counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1).fillna(0) * 100
    order = rel.sum(axis=1).sort_values(ascending=False).index.tolist()
    top = order[:top_n]
    shown = rel.loc[top].copy()
    rest = rel.drop(index=top, errors="ignore").sum(axis=0)
    if rest.sum() > 0:
        shown.loc["Other"] = rest
    return shown


def stacked_bar(ax, counts: pd.DataFrame, top_n: int, title: str, ylabel: str = "Relative abundance (%)"):
    rel = relative_top_table(counts, top_n)
    sites = ["Lung", "Trachea"]
    bottom = np.zeros(len(sites))
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(rel), 1)))
    for i, (taxon, row) in enumerate(rel.iterrows()):
        vals = row.reindex(sites).fillna(0).values
        ax.bar(sites, vals, bottom=bottom, label=taxon, color=colors[i])
        bottom += vals
    ax.set_ylim(0, 100)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, frameon=False)
    ax.grid(axis="y", alpha=0.15)


def savefig_all(fig, stem: str, dpi: int = 600):
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(FIGURES_DIR / f"{stem}.{ext}", dpi=dpi, bbox_inches="tight")


def draw_venn(ax, lung_set: set[str], trachea_set: set[str]):
    lung_only = len(lung_set - trachea_set)
    shared = len(lung_set & trachea_set)
    trachea_only = len(trachea_set - lung_set)
    if venn2 is not None:
        venn2(subsets=(lung_only, trachea_only, shared), set_labels=("Lung", "Trachea"), ax=ax)
    else:
        from matplotlib.patches import Circle
        ax.add_patch(Circle((0.42, 0.5), 0.28, alpha=0.35))
        ax.add_patch(Circle((0.58, 0.5), 0.28, alpha=0.35))
        ax.text(0.30, 0.50, str(lung_only), ha="center", va="center", fontsize=14)
        ax.text(0.50, 0.50, str(shared), ha="center", va="center", fontsize=14, fontweight="bold")
        ax.text(0.70, 0.50, str(trachea_only), ha="center", va="center", fontsize=14)
        ax.text(0.30, 0.18, "Lung", ha="center")
        ax.text(0.70, 0.18, "Trachea", ha="center")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
    ax.set_title("Genus presence", fontweight="bold")


def heatmap(ax, genus_counts: pd.DataFrame):
    lung_set = set(genus_counts.index[genus_counts["Lung"] > 0])
    trachea_set = set(genus_counts.index[genus_counts["Trachea"] > 0])
    total = genus_counts.sum(axis=1)
    groups = [
        sorted(lung_set - trachea_set, key=lambda g: total[g], reverse=True),
        sorted(lung_set & trachea_set, key=lambda g: total[g], reverse=True),
        sorted(trachea_set - lung_set, key=lambda g: total[g], reverse=True),
    ]
    order = [g for group in groups for g in group]
    mat_counts = genus_counts.loc[order, ["Lung", "Trachea"]].T
    mat = np.where(mat_counts.values > 0, np.log10(mat_counts.values), -1.0)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under("#f3e7c0")
    im = ax.imshow(mat, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0)
    ax.set_yticks(range(mat_counts.shape[0]))
    ax.set_yticklabels(mat_counts.index)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, ha="center", fontsize=7)
    ax.set_title("Genus-level counts by body site", fontweight="bold")
    ax.set_xlabel("Genus")
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("pseudo-log10 count; beige = 0")
    return order


def plot_shannon(metadata: pd.DataFrame):
    if not SHANNON.exists():
        print(f"[WARN] Shannon file not found: {SHANNON}")
        return
    sh = pd.read_csv(SHANNON, sep="\t", dtype={0: str})
    sh = sh.rename(columns={sh.columns[0]: "sample-id", sh.columns[1]: "Shannon"})
    sh["Shannon"] = pd.to_numeric(sh["Shannon"], errors="coerce")
    data = sh.merge(metadata[["sample-id", "body_site"]], on="sample-id", how="inner").dropna()
    groups = [data.loc[data["body_site"] == site, "Shannon"].values for site in ["Lung", "Trachea"]]
    if all(len(g) > 0 for g in groups):
        stat, pval = kruskal(*groups)
    else:
        stat, pval = math.nan, math.nan
    pd.DataFrame([{"test": "Kruskal-Wallis", "H": stat, "p_value": pval}]).to_csv(
        EXPORTED_DIR / "shannon_kruskal_wallis.tsv", sep="\t", index=False
    )
    fig, ax = plt.subplots(figsize=(4.2, 4.5))
    values = [data.loc[data["body_site"] == site, "Shannon"].values for site in ["Lung", "Trachea"]]
    ax.boxplot(values, labels=["Lung", "Trachea"], showmeans=True)
    for i, vals in enumerate(values, start=1):
        ax.scatter(np.repeat(i, len(vals)), vals, s=30, zorder=3)
    ax.set_ylabel("Shannon diversity (H)")
    ax.set_title(f"Rarefied Shannon diversity\nKruskal-Wallis H={stat:.3g}, p={pval:.3g}")
    ax.grid(axis="y", alpha=0.2)
    savefig_all(fig, "shannon_alpha_diversity")
    plt.close(fig)


def main():
    metadata = read_metadata(METADATA)
    class_counts = group_by_site(read_biom_tsv(TAXA_L3), metadata)
    genus_counts = group_by_site(read_biom_tsv(TAXA_L6), metadata)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    stacked_bar(ax, class_counts, top_n=8, title="A. Bacterial classes")
    savefig_all(fig, "figure3A_class_barplot")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    stacked_bar(ax, genus_counts, top_n=14, title="B. Bacterial genera")
    savefig_all(fig, "figure3B_genus_barplot")
    plt.close(fig)

    lung_set = set(genus_counts.index[genus_counts["Lung"] > 0])
    trachea_set = set(genus_counts.index[genus_counts["Trachea"] > 0])
    fig, ax = plt.subplots(figsize=(4.5, 4.0))
    draw_venn(ax, lung_set, trachea_set)
    savefig_all(fig, "figure3C_venn")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(max(8, len(genus_counts) * 0.32), 3.8))
    order = heatmap(ax, genus_counts)
    savefig_all(fig, "figure3D_genus_heatmap")
    plt.close(fig)

    genus_counts.loc[order].to_csv(EXPORTED_DIR / "genus_counts_lung_trachea_ordered.tsv", sep="\t")

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.9, 1.25], wspace=0.55, hspace=0.75)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    ax_d = fig.add_subplot(gs[2, :])
    stacked_bar(ax_a, class_counts, top_n=8, title="A. Bacterial classes")
    stacked_bar(ax_b, genus_counts, top_n=14, title="B. Bacterial genera")
    draw_venn(ax_c, lung_set, trachea_set)
    heatmap(ax_d, genus_counts)
    fig.suptitle("Composition of broiler respiratory microbiomes", fontweight="bold", y=0.995)
    savefig_all(fig, "Figure_3_panel")
    plt.close(fig)

    plot_shannon(metadata)
    print(f"Figures saved in {FIGURES_DIR}")


if __name__ == "__main__":
    main()
