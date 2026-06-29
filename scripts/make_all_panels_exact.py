#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Figura A4 com 4 painéis:
  A) Classe (top-5 + Others + Unclassified [= Unassigned + sem rótulo no nível])
  B) Gênero (top-12 + Others + Unclassified [= Unassigned + sem rótulo no nível])
  C) Venn (Lung × Trachea) – à esquerda, abaixo de A
  D) Heatmap por gênero (mesmo conjunto de C) em escala pseudo-log10:
     - zeros mapeados para -1 (exibidos como 10^-1 na colorbar, cor bege)
     - contagens > 0 mapeadas para log10(count)
     - colorbar com ticks em 10^-1, 10^0, 10^1, … até o expoente máximo

Saídas em --outdir:
  fig_ABCD_A4.(png|svg), venn.(png|svg), heatmap_genus.(png|svg)
  Tabelas: debug_class_relative.csv, debug_genus_relative.csv,
           genus_presence_counts.tsv, heatmap_matrix_by_genus.tsv,
           genus_table_long.tsv, diffab_features.tsv, diffab_genus.tsv,
           genus_{Lung|Trachea}_only.txt, genus_both.txt
"""

import os, re, math, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter, MultipleLocator, FixedLocator, FixedFormatter
from matplotlib.patches import Circle
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu, norm as stats_norm
from qiime2 import Artifact
from biom import Table

# ---------- Paleta ----------
COLOR_UNASSIGNED = "#E0E0E0"   # cor usada para o grupo Unclassified (fundido)
COLOR_OTHERS     = "#9E9E9E"
FALLBACK = [
    "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
    "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf",
    "#393b79","#637939","#8c6d31","#843c39","#7b4173",
    "#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5"
]

VENN_RADIUS = 1.8
VENN_GAP    = 2.2
VENN_FONT   = 14

# ---------- Args ----------
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="FeatureTable[Frequency] .qza (NÃO agrupada por site)")
    ap.add_argument("--taxonomy", required=True, help="FeatureData[Taxonomy] .qza")
    ap.add_argument("--metadata", required=True, help="sample-metadata.tsv")
    ap.add_argument("--group-col", default="Site", help="Coluna de agrupamento (ex.: Site)")
    ap.add_argument("--group-a", default="Lung", help="Nome do grupo A")
    ap.add_argument("--group-b", default="Trachea", help="Nome do grupo B")
    ap.add_argument("--topA", type=int, default=5,  help="Top-N Classe (A)")
    ap.add_argument("--topB", type=int, default=12, help="Top-N Gênero (B)")
    ap.add_argument("--alpha", type=float, default=0.05, help="FDR (BH) (tabelas)")
    ap.add_argument("--min-count", type=float, default=1.0, help="Presença: contagem mínima em ≥ min-samples")
    ap.add_argument("--min-samples", type=int, default=1, help="Presença: nº mínimo de amostras no grupo")
    ap.add_argument("--pseudocount", type=float, default=0.5, help="Pseudocontagem (CLR / proporções)")
    ap.add_argument("--heatmap-top", type=int, default=60, help="Máximo de gêneros exibidos no Venn/Heatmap")
    ap.add_argument("--outdir", required=True)
    return ap.parse_args()

# ---------- Utils ----------
def biom_to_df(tbl: Table) -> pd.DataFrame:
    return pd.DataFrame(tbl.matrix_data.toarray(),
                        index=tbl.ids("observation"),
                        columns=tbl.ids("sample"))

def load_metadata(path):
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    if 'sample-id' not in df.columns:
        df = df.rename(columns={df.columns[0]: 'sample-id'})
    return df

def truncate_to_level(taxon: str, level: str) -> str:
    """Trunca a linhagem ao nível (phylum=2, class=3, genus=6).
       Se 'Unassigned', gera 'Unassigned;__;...'
    """
    if not isinstance(taxon, str) or taxon.startswith("Unassigned"):
        depth = {"phylum": 2, "class": 3, "genus": 6}[level]
        return ";".join(["Unassigned"] + ["__"]*(depth-1))
    parts = [p.strip() for p in taxon.split(";")]
    need  = {"phylum": 2, "class": 3, "genus": 6}[level]
    if len(parts) < need:
        parts += ["__"]*(need-len(parts))
    return ";".join(parts[:need])

def last_label(lineage: str) -> str:
    parts = [p.strip() for p in str(lineage).split(";")]
    return parts[-1] if parts else "__"

def is_blank_label(label: str) -> bool:
    if label is None: return True
    s = str(label).strip()
    if s in ("", "__"): return True
    if re.fullmatch(r"[a-zA-Z]__", s): return True  # g__, f__, etc.
    return False

def is_unassigned_lineage(lineage: str) -> bool:
    return str(lineage).strip().startswith("Unassigned")

def is_unclassified_label(label: str) -> bool:
    """Tudo que deve entrar no grupo Unclassified (no nível visual)."""
    return is_blank_label(label) or str(label).strip() == "Unassigned"

def strip_prefix(label: str) -> str:
    return re.sub(r"^[a-zA-Z]__", "", str(label)).strip()

def pretty_label_class(lab: str) -> str:
    if lab == "Unclassified": return "Unclassified"
    if lab == "Others": return "Others"
    return strip_prefix(lab)

def pretty_label_genus(lab: str) -> str:
    if lab == "Unclassified": return "Unclassified"
    if lab == "Others": return "Others"
    return strip_prefix(lab)

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    m = len(p); order = np.argsort(p)
    q = np.empty(m, dtype=float); prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        i = m - rank + 1
        val = min(prev, p[idx] * m / i)
        q[idx] = val; prev = val
    return np.clip(q, 0, 1)

def clr(x, pseudocount=0.5):
    x = np.asarray(x, dtype=float) + pseudocount
    lx = np.log(x)
    return lx - lx.mean()

# ---------- Taxonomy aggregation ----------
def aggregate_level(counts_df, taxonomy_series, level):
    """Retorna tabela agregada por label visual do nível:
       - Unassigned + rank vazio => Unclassified
       - demais: último rótulo do nível (ex.: c__Gammaproteobacteria)
    """
    mapping = {}
    for fid in counts_df.index:
        lin = taxonomy_series.get(fid, "Unassigned")
        lin_trunc = truncate_to_level(lin, level)
        lab = last_label(lin_trunc)
        if is_unassigned_lineage(lin_trunc) or is_unclassified_label(lab):
            key = "Unclassified"
        else:
            key = lab
        mapping[fid] = key
    tmp = counts_df.copy()
    tmp["__tax__"] = [mapping[i] for i in tmp.index]
    return tmp.groupby("__tax__").sum()

def choose_top_with_unclassified_and_others(rel_by_group, topn):
    """rel_by_group: rows=labels, cols=groups, valores %.
       Mantém Unclassified sempre separado, escolhe topn entre os demais por média,
       soma o resto como Others.
    """
    idx = list(rel_by_group.index)
    has_un = "Unclassified" in idx
    non_un = [i for i in idx if i != "Unclassified"]
    mean_abund = rel_by_group.loc[non_un].mean(axis=1).sort_values(ascending=False) if non_un else pd.Series(dtype=float)
    top = list(mean_abund.head(topn).index)
    rest = [i for i in non_un if i not in top]
    rows = []
    if has_un:
        rows.append("Unclassified")
    rows += top
    out = rel_by_group.loc[rows].copy() if rows else pd.DataFrame(columns=rel_by_group.columns)
    if rest:
        out.loc["Others"] = rel_by_group.loc[rest].sum(axis=0)
    return out

# ---------- Plots ----------
def _draw_stack(ax, df_percent, labels_order, pretty_func, is_genus=False, legend_side="right", add_ylabel=True, ticks_left=True):
    groups = list(df_percent.columns)
    bottoms = np.zeros(len(groups))
    color_map = {}
    # cor fixa para Unclassified e Others; demais paleta
    k = 0
    for lab in labels_order:
        if lab == "Unclassified": color_map[lab] = COLOR_UNASSIGNED
        elif lab == "Others":     color_map[lab] = COLOR_OTHERS
        else:
            color_map[lab] = FALLBACK[k % len(FALLBACK)]; k += 1
    for lab in labels_order:
        vals = df_percent.loc[lab, groups].values.astype(float)
        ax.bar(groups, vals, bottom=bottoms, color=color_map[lab], edgecolor="white", linewidth=0.7, width=0.66, label=pretty_func(lab))
        bottoms += vals
    ax.set_ylim(0, 100)
    if add_ylabel:
        ax.set_ylabel("Relative Abundance (%)", fontsize=10)
    ax.yaxis.set_major_formatter(PercentFormatter(100))
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=9, labelleft=ticks_left)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y', linestyle=':', alpha=0.28)
    ax.set_axisbelow(True)
    # legenda
    handles, labels = ax.get_legend_handles_labels()
    if legend_side == "left":
        ax.legend(handles, labels, loc="center right", bbox_to_anchor=(-0.18, 0.50),
                  frameon=False, fontsize=7.2, title="Class" if not is_genus else "Genus",
                  title_fontsize=8.5, borderaxespad=0.0)
    else:
        ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.50),
                  frameon=False, fontsize=6.8 if is_genus else 7.2,
                  title="Genus" if is_genus else "Class", title_fontsize=8.5,
                  borderaxespad=0.0)

# ---------- Venn sem matplotlib_venn ----------
def draw_venn_on_axis(ax, nAonly, nBoth, nBonly, labels=("Lung","Trachea")):
    ax.set_aspect('equal'); ax.axis('off')
    r = VENN_RADIUS
    c1 = np.array([0.0, 0.0]); c2 = np.array([VENN_GAP, 0.0])
    ax.add_patch(Circle(c1, r, facecolor="#9BC8E6", edgecolor="black", lw=1.0, alpha=0.75))
    ax.add_patch(Circle(c2, r, facecolor="#F4A3A3", edgecolor="black", lw=1.0, alpha=0.75))
    ax.text(c1[0]-0.72, 0, str(nAonly), ha='center', va='center', fontsize=VENN_FONT, fontweight='bold')
    ax.text((c1[0]+c2[0])/2, 0, str(nBoth), ha='center', va='center', fontsize=VENN_FONT, fontweight='bold')
    ax.text(c2[0]+0.72, 0, str(nBonly), ha='center', va='center', fontsize=VENN_FONT, fontweight='bold')
    ax.text(c1[0], -r-0.50, labels[0], ha='center', va='top', fontsize=11)
    ax.text(c2[0], -r-0.50, labels[1], ha='center', va='top', fontsize=11)
    ax.set_xlim(-r-0.65, VENN_GAP+r+0.65)
    ax.set_ylim(-r-0.95, r+0.45)

# ---------- Heatmap ----------
def make_heat_cmap():
    base = cm.get_cmap("YlOrRd", 256)
    colors = base(np.linspace(0, 1, 256))
    cmap = LinearSegmentedColormap.from_list("beige_YlOrRd", colors)
    cmap.set_under("#F3E6C8")  # zeros
    return cmap

def draw_heatmap_on_axis(ax, counts_by_genus_site, labels_x, gA, gB, splitA, splitBoth):
    # matriz 2 x N com pseudo-log10: zero => -1; >0 => log10(count)
    mat_counts = counts_by_genus_site.loc[labels_x, [gA, gB]].T.astype(float)
    mat = mat_counts.values.copy()
    mat_log = np.where(mat > 0, np.log10(mat), -1.0)
    vmax = max(0.0, np.nanmax(mat_log))
    im = ax.imshow(mat_log, aspect='auto', interpolation='nearest', cmap=make_heat_cmap(), vmin=0, vmax=vmax)
    ax.set_yticks([0,1]); ax.set_yticklabels([gA, gB], fontsize=10)
    ax.set_xticks(np.arange(len(labels_x)))
    ax.set_xticklabels([strip_prefix(x) for x in labels_x], rotation=60, ha='right', fontsize=7)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(axis='y', length=0)
    # linhas separadoras Lung-only | Both | Trachea-only
    if splitA > 0:
        ax.axvline(splitA-0.5, color='black', lw=0.8)
    if splitA + splitBoth > 0:
        ax.axvline(splitA+splitBoth-0.5, color='black', lw=0.8)
    # colorbar com ticks pseudo-log
    cbar = plt.colorbar(im, ax=ax, fraction=0.018, pad=0.008)
    max_exp = int(math.ceil(max(0, vmax)))
    ticks = [-1.0] + list(range(0, max_exp+1))
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([r"$10^{-1}$ (0)"] + [rf"$10^{i}$" for i in range(0, max_exp+1)])
    cbar.set_label("Pseudo-log10(count)", fontsize=9)
    ax.spines[['top','right','left','bottom']].set_visible(False)

# ---------- Differential exploratory tables ----------
def diff_tables(counts_feature, taxonomy_series, md, group_col, gA, gB, pseudocount=0.5, alpha=0.05):
    md2 = md.set_index('sample-id')
    samplesA = [s for s in counts_feature.columns if s in md2.index and md2.loc[s, group_col] == gA]
    samplesB = [s for s in counts_feature.columns if s in md2.index and md2.loc[s, group_col] == gB]
    rows = []
    # feature-level proportions + CLR, Mann-Whitney simples (n pequeno -> exploratório)
    lib = counts_feature.sum(axis=0).replace(0, np.nan)
    prop = counts_feature.divide(lib, axis=1).fillna(0.0)
    for fid in counts_feature.index:
        a = prop.loc[fid, samplesA].values
        b = prop.loc[fid, samplesB].values
        try:
            p = mannwhitneyu(a, b, alternative='two-sided').pvalue
        except Exception:
            p = 1.0
        meanA = a.mean() if len(a) else 0
        meanB = b.mean() if len(b) else 0
        lfc = np.log2((meanA+pseudocount/lib[samplesA].mean())/(meanB+pseudocount/lib[samplesB].mean())) if len(samplesA) and len(samplesB) else np.nan
        rows.append([fid, taxonomy_series.get(fid, "Unassigned"), meanA, meanB, lfc, p])
    df = pd.DataFrame(rows, columns=['feature_id','taxonomy','mean_rel_'+gA,'mean_rel_'+gB,'log2FC_'+gA+'_vs_'+gB,'p'])
    df['q_BH'] = bh_fdr(df['p'].values) if len(df) else []
    return df.sort_values(['q_BH','p'])

def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    table_art = Artifact.load(a.table)
    tax_art   = Artifact.load(a.taxonomy)
    md = load_metadata(a.metadata)
    # Accept both Site and body_site metadata conventions
    if a.group_col not in md.columns and a.group_col == "Site" and "body_site" in md.columns:
        a.group_col = "body_site"
    if a.group_col not in md.columns:
        raise SystemExit(f"Coluna '{a.group_col}' não encontrada no metadata. Colunas: {list(md.columns)}")

    counts = biom_to_df(table_art.view(Table))
    tax = tax_art.view(pd.DataFrame)
    if 'Taxon' not in tax.columns:
        raise SystemExit("taxonomy.qza não tem coluna 'Taxon'")
    taxonomy = tax['Taxon']

    # Alinhar amostras
    md_samples = set(md['sample-id'])
    samples = [s for s in counts.columns if s in md_samples]
    counts = counts[samples]
    md = md[md['sample-id'].isin(samples)].copy()

    gA, gB = a.group_a, a.group_b
    md2 = md.set_index('sample-id')
    samplesA = [s for s in samples if md2.loc[s, a.group_col] == gA]
    samplesB = [s for s in samples if md2.loc[s, a.group_col] == gB]
    if len(samplesA) == 0 or len(samplesB) == 0:
        raise SystemExit(f"Amostras insuficientes nos grupos: {gA}={len(samplesA)}, {gB}={len(samplesB)}")

    # A/B: agregar por classe e gênero
    class_counts = aggregate_level(counts, taxonomy, "class")
    genus_counts = aggregate_level(counts, taxonomy, "genus")

    def group_sum(tbl):
        return pd.DataFrame({gA: tbl[samplesA].sum(axis=1),
                             gB: tbl[samplesB].sum(axis=1)})
    class_site = group_sum(class_counts)
    genus_site = group_sum(genus_counts)

    class_rel = class_site.divide(class_site.sum(axis=0), axis=1)*100.0
    genus_rel = genus_site.divide(genus_site.sum(axis=0), axis=1)*100.0
    A_plot = choose_top_with_unclassified_and_others(class_rel, a.topA)
    B_plot = choose_top_with_unclassified_and_others(genus_rel, a.topB)
    labels_A = list(A_plot.index)
    labels_B = list(B_plot.index)

    # C/D presença por gênero (Unclassified tratado como um gênero visual)
    genus_sample = genus_counts.copy()
    presenceA = (genus_sample[samplesA] >= a.min_count).sum(axis=1) >= a.min_samples
    presenceB = (genus_sample[samplesB] >= a.min_count).sum(axis=1) >= a.min_samples
    cats = pd.Series(index=genus_sample.index, dtype=str)
    cats[presenceA & ~presenceB] = f"{gA}-only"
    cats[presenceA & presenceB]  = "Both"
    cats[~presenceA & presenceB] = f"{gB}-only"
    cats = cats.dropna()

    heat_counts_full = genus_site.loc[cats.index]
    # ordenação: Lung-only por contagem em A, Both por soma, Trachea-only por contagem em B
    Aonly = list(cats[cats.eq(f"{gA}-only")].index)
    Both  = list(cats[cats.eq("Both")].index)
    Bonly = list(cats[cats.eq(f"{gB}-only")].index)
    Aonly = list(heat_counts_full.loc[Aonly, gA].sort_values(ascending=False).index)
    Both  = list(heat_counts_full.loc[Both].sum(axis=1).sort_values(ascending=False).index)
    Bonly = list(heat_counts_full.loc[Bonly, gB].sort_values(ascending=False).index)
    labels_x = Aonly + Both + Bonly
    # limitar heatmap se necessário, preservando proporção por blocos simples
    if len(labels_x) > a.heatmap_top:
        labels_x = labels_x[:a.heatmap_top]
    cats = cats.loc[labels_x]
    heat_counts = heat_counts_full.loc[labels_x]
    nAonly = int((cats == f"{gA}-only").sum())
    nBoth  = int((cats == "Both").sum())
    nBonly = int((cats == f"{gB}-only").sum())

    # tabelas de apoio
    presence_table = pd.DataFrame({"genus": labels_x, "category": cats.values,
                                   gA: heat_counts[gA].values, gB: heat_counts[gB].values})
    presence_table.to_csv(os.path.join(a.outdir, "genus_presence_counts.tsv"), sep="\t", index=False)
    heat_counts.to_csv(os.path.join(a.outdir, "heatmap_matrix_by_genus.tsv"), sep="\t")
    genus_site.reset_index().rename(columns={'__tax__':'genus'}).to_csv(os.path.join(a.outdir, "genus_table_long.tsv"), sep="\t", index=False)

    # diff exploratory: feature e gênero
    df_feat = diff_tables(counts, taxonomy, md, a.group_col, gA, gB, a.pseudocount, a.alpha)
    df_feat.to_csv(os.path.join(a.outdir, "diffab_features.tsv"), sep="\t", index=False)
    # gênero: reutiliza genus_sample como "features"
    fake_tax = pd.Series(genus_sample.index, index=genus_sample.index)
    df_gen = diff_tables(genus_sample, fake_tax, md, a.group_col, gA, gB, a.pseudocount, a.alpha)
    df_gen = df_gen.rename(columns={'feature_id':'genus','taxonomy':'genus_label'})
    df_gen.to_csv(os.path.join(a.outdir, "diffab_genus.tsv"), sep="\t", index=False)

    # -------- Layout A4 landscape --------
    fig_w_in, fig_h_in = 11.69, 8.27
    left_m, right_m = 1.15, 1.15
    top_m = 0.45
    gap_in = 1.35
    wA = (fig_w_in - left_m - right_m - gap_in) / 2
    wB = wA
    hA = 2.35
    hC = 1.90
    hD = 1.65
    gv1, gv2 = 0.28, 0.35

    leftA = left_m / fig_w_in
    widthA  = wA / fig_w_in
    leftB   = (left_m + wA + gap_in) / fig_w_in
    widthB  = wB / fig_w_in
    bottomA = (fig_h_in - top_m - hA) / fig_h_in
    bottomB = bottomA
    heightAB = hA / fig_h_in

    leftC  = leftA; widthC = widthA
    bottomC = (fig_h_in - top_m - hA - gv1 - hC) / fig_h_in
    heightC = hC / fig_h_in

    leftD  = left_m / fig_w_in
    widthD = (fig_w_in - left_m - right_m) / fig_w_in
    bottomD = (fig_h_in - top_m - hA - gv1 - hC - gv2 - hD) / fig_h_in
    heightD = hD / fig_h_in

    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=300)

    # Painel A
    axA = fig.add_axes([leftA, bottomA, widthA, heightAB])
    _draw_stack(axA, A_plot, labels_A, pretty_label_class,
                is_genus=False, legend_side="left", add_ylabel=True, ticks_left=True)
    axA.text(-0.12, 1.12, "A", transform=axA.transAxes, fontsize=16, fontweight="bold", va="top")

    # Painel B
    axB = fig.add_axes([leftB, bottomB, widthB, heightAB])
    _draw_stack(axB, B_plot, labels_B, pretty_label_genus,
                is_genus=True, legend_side="right", add_ylabel=False, ticks_left=True)
    axB.text(-0.12, 1.12, "B", transform=axB.transAxes, fontsize=16, fontweight="bold", va="top")

    # Painel C
    axC = fig.add_axes([leftC, bottomC, widthC, heightC])
    draw_venn_on_axis(axC, nAonly, nBoth, nBonly, labels=(gA, gB))
    axC.text(-0.08, 1.08, "C", transform=axC.transAxes, fontsize=16, fontweight="bold", va="top")

    # Painel D — heatmap pseudo-log10
    axD = fig.add_axes([leftD, bottomD, widthD, heightD])
    draw_heatmap_on_axis(axD, heat_counts, labels_x, gA, gB, nAonly, nBoth)
    axD.text(-0.05, 1.06, "D", transform=axD.transAxes, fontsize=16, fontweight="bold", va="top")

    fig.savefig(os.path.join(a.outdir, "fig_ABCD_A4.png"), dpi=1200, bbox_inches="tight")
    fig.savefig(os.path.join(a.outdir, "fig_ABCD_A4.svg"), bbox_inches="tight")
    plt.close(fig)

    # Versões individuais (opcional)
    figC, axC2 = plt.subplots(figsize=(8.0, 6.0), dpi=300)
    draw_venn_on_axis(axC2, nAonly, nBoth, nBonly, labels=(gA, gB))
    axC2.text(-0.08, 1.08, "C", transform=axC2.transAxes, fontsize=16, fontweight="bold", va="top")
    figC.tight_layout(); figC.savefig(os.path.join(a.outdir, "venn.png"), dpi=1200, bbox_inches="tight")
    figC.savefig(os.path.join(a.outdir, "venn.svg"), bbox_inches="tight"); plt.close(figC)

    ncols = len(labels_x)
    fig_w = max(9.0, 0.35 * ncols + 6.0)
    fig_h = 2.4
    figD, axD2 = plt.subplots(1, 1, figsize=(fig_w, fig_h), dpi=300)
    draw_heatmap_on_axis(axD2, heat_counts, labels_x, gA, gB, nAonly, nBoth)
    axD2.text(-0.05, 1.06, "D", transform=axD2.transAxes, fontsize=16, fontweight="bold", va="top")
    figD.tight_layout()
    figD.savefig(os.path.join(a.outdir, "heatmap_genus.png"), dpi=1200, bbox_inches="tight")
    figD.savefig(os.path.join(a.outdir, "heatmap_genus.svg"), bbox_inches="tight")
    plt.close(figD)

    # Listas de gêneros por categoria
    with open(os.path.join(a.outdir, f"genus_{gA}_only.txt"), "w") as f:
        for g in labels_x[cats.eq(f"{gA}-only")]: f.write(g + "\n")
    with open(os.path.join(a.outdir, "genus_both.txt"), "w") as f:
        for g in labels_x[cats.eq("Both")]: f.write(g + "\n")
    with open(os.path.join(a.outdir, f"genus_{gB}_only.txt"), "w") as f:
        for g in labels_x[cats.eq(f"{gB}-only")]: f.write(g + "\n")

    # depuração A/B
    A_plot.to_csv(os.path.join(a.outdir, "debug_class_relative.csv"))
    B_plot.to_csv(os.path.join(a.outdir, "debug_genus_relative.csv"))

    print(f"[INFO] Heatmap (pseudo-log10): max_count={heat_counts.max():.3f} -> vmax_log={np.log10(max(1.0, heat_counts.max())):.3f}")
    print("[OK] Figura completa: fig_ABCD_A4.(png|svg)")
    print("[OK] Tabelas salvas em:", os.path.abspath(a.outdir))

if __name__ == "__main__":
    main()
