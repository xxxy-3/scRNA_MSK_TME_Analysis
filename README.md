# Single-Cell Transcriptomics of the Lung Tumour Microenvironment

**GSE123902 · primary lung adenocarcinoma vs. non-tumour lung · 41,384 cells · 12 samples**

A reproducible Scanpy pipeline that goes from raw per-sample `.h5ad` files to annotated cell types, differential expression, cell–cell communication and pseudotime — built as an end-to-end exercise in single-cell analysis rather than a narrow re-implementation of one figure.

---

## Overview

The tumour microenvironment (TME) of lung adenocarcinoma differs from adjacent non-tumour lung not only in its cell-type composition but also in the transcriptional state of the cells that are present. This repository asks four questions of a public 10x Genomics dataset:

1. **Composition** — which cell types make up tumour vs. non-tumour tissue, and how variable is that between patients?
2. **Transcriptional shift** — which genes distinguish tumour from non-tumour cells once samples are integrated?
3. **Communication** — which ligand–receptor programmes connect the cell types in this ecosystem?
4. **Trajectory** — do cells from the two conditions occupy different positions along a differentiation axis?

Everything below is produced by one notebook (`scRNA_pipeline.ipynb`) plus the figures in `figures/`.

---

## Dataset

| Item | Value |
|---|---|
| Source | [GEO: GSE123902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123902) |
| Tissue | Primary lung adenocarcinoma and non-tumour-involved lung |
| Samples used | **8 primary tumour** (`LX653/LX661/LX675/LX676/LX679/LX680/LX682/LX684_T`) + **4 non-tumour** (`LX675/LX682/LX684/LX685_N`) |
| Platform | 10x Genomics Chromium, 3′ scRNA-seq |
| Cells profiled in the series | 41,384 (the series additionally contains brain, bone and adrenal metastases, which are **excluded** here) |
| Input format | One `.h5ad` file per sample (not committed — see `.gitignore`) |

Only samples whose names contain `PRIMARY_TUMOUR` or `NORMAL` are read; metastatic samples are filtered out at load time so that the comparison stays tumour-vs-normal within the lung.

---

## Results

### 1. Cell-type composition differs between tumour and non-tumour tissue

Cells were annotated with CellTypist (`Immune_All_Low`, majority voting) and collapsed into five major immune compartments plus the non-immune labels returned by the model (epithelial, endothelial, mast cells).

![Cell type fractions per sample](figures/vertical_barplot_fixed.png)

T cells are the largest compartment in almost every sample. The clearest difference is epithelial content: `LX680_T`, `LX684_T`, `LX679_T` and `LX653_T` carry a substantial epithelial fraction (>20%), whereas all four non-tumour samples are dominated by macrophages/monocytes and T cells with little epithelium. Composition is markedly patient-specific — `LX676_T` is almost entirely immune — which is itself the main reason the differential-expression results below should be read with care.

### 2. Tumour cells show a chaperone-high, immune-effector-low transcriptional programme

Wilcoxon rank-sum test on all genes, `Primary Tumor` vs `Non-tumour` (reference), BH-corrected.

![Volcano plot](figures/volcano_plot.png)

Thresholds: |log₂FC| > 1 and adjusted *P* < 0.05. The distribution is strongly asymmetric — many more genes reach significance on the down-regulated side at small effect sizes, while a smaller set of up-regulated genes reaches very large fold changes (log₂FC up to ~8).

![Top differentially expressed genes](figures/DEG_heatmap_NCS_final.png)

Top 25 up- and 25 down-regulated genes (by adjusted *P*), shown for 150 randomly sampled cells per group, row-wise Z-scored and clipped at ±2. Two coherent programmes stand out:

- **Up in tumour** — a heat-shock / chaperone module (`HSPA1A`, `HSP90AA1`, `HSPH1`, `HSPD1`, `HSPE1`, `HSPA8`, `HSPB1`, `DNAJB1`, `DNAJA1`), together with ribosomal proteins (`RPLP1`, `RPL28`, `RPL36`) and immediate-early / stress-response genes (`JUN`, `DUSP4`, `NR4A1`, `BTG1`).
- **Down in tumour** — immune-effector and myeloid adaptor transcripts (`NKG7`, `CCL4`, `FCER1G`, `TYROBP`, `ITGB2`, `CD52`, `HLA-C`), plus alveolar/lipid-associated genes (`SCGB1A1`, `FABP4`, `FTL`).

The down-regulated side is consistent with the composition shift in Figure 1 (fewer NK/T and myeloid cells per cell sampled). The up-regulated chaperone module is *not* composition-driven in the same way and is discussed as a caveat in [Limitations](#limitations).

### 3. The TME is densely connected, with macrophages and epithelium as hubs

Consensus ligand–receptor inference with LIANA across the nine annotated populations, restricted to high-confidence interactions (`magnitude_rank` ≤ 0.01).

![Cell-cell communication](figures/cell_communication_circle.png)

Interactions are near-saturating in this dataset: epithelial cells, macrophages/monocytes, fibroblasts-like stromal labels and endothelial cells all act as both senders and receivers. Edge density of this kind should be read as "these populations are transcriptionally competent to signal to each other", not as evidence that a given interaction occurs *in situ* — LIANA ranks expression of known pairs, it does not measure contact.

### 4. Pseudotime separates the two conditions only partially

Diffusion pseudotime (DPT) on the Harmony-corrected PCA space, with the root set to a cell in Leiden cluster 0.

Pseudotime distributions for the two groups overlap substantially, i.e. the dominant axis of variation in this dataset is **cell identity, not disease status**. A meaningful trajectory analysis here would require restricting to one lineage (e.g. the monocyte→macrophage axis) and choosing the root from marker expression rather than from an arbitrary cluster — this is listed under [What I would do next](#what-i-would-do-next).

---

## Pipeline

| Step | Method | Key parameters |
|---|---|---|
| Load & merge | `anndata.concat` over per-sample files | `index_unique="-"`, batch key = sample name |
| Cell QC | `filter_cells` / `filter_genes` | ≥200 genes per cell, gene in ≥3 cells |
| Mitochondrial QC | `calculate_qc_metrics` | retain `pct_counts_mt < 15` |
| Normalisation | `normalize_total` + `log1p` | target sum 1e4; raw counts kept in `.raw` |
| HVG selection | `highly_variable_genes` | `min_mean=0.0125`, `max_mean=3`, `min_disp=0.5`, batch-aware |
| Scaling / PCA | `scale` + `pca` | `max_value=10`, `svd_solver="arpack"` |
| Batch correction | **Harmony** (`harmonypy`, called directly) | on `batch`; corrected embedding in `obsm["X_pca_harmony"]` |
| Neighbours / UMAP | `pp.neighbors` + `tl.umap` | `n_neighbors=15`, `n_pcs=30` |
| Clustering | Leiden | `resolution=0.6`, `flavor="igraph"`, `n_iterations=2` |
| Annotation | **CellTypist** `Immune_All_Low` + majority voting | manual collapse to 5 major immune types |
| Composition | fraction per sample, stacked bar | tumour and normal ordered separately |
| DE | `rank_genes_groups`, Wilcoxon | reference = `Non-tumour`; exported to CSV |
| Visualisation | volcano + clustered heatmap | \|log₂FC\|>1, adj. *P*<0.05; Z-score clipped ±2 |
| Communication | **LIANA** consensus + circle plot | `magnitude_rank ≤ 0.01` |
| Trajectory | `diffmap` + `dpt` | root = first cell of Leiden cluster 0 |

`harmonypy` is called directly rather than through `sc.external.pp.harmony_integrate`, because the Scanpy wrapper returns an AnnData view whose sparse matrix triggers `getnnz` errors on some versions; the corrected matrix is transposed if necessary and cast to `float32`.

---

## Repository structure

```
scRNA_MSK_TME_Analysis/
├── scRNA_pipeline.ipynb          # the full analysis, cell by cell
├── requirements.txt
├── figures/
│   ├── vertical_barplot_fixed.pdf/.png   # composition, per sample
│   ├── cell_type_fractions.pdf/.png      # composition, grouped bars
│   ├── volcano_plot.pdf/.png             # DE, tumour vs non-tumour
│   ├── DEG_heatmap_NCS_final.pdf/.png    # top 50 DEGs
│   └── cell_communication_circle.pdf/.png
└── README.md
```

PDFs are the vector originals; PNGs at 200 dpi exist so that the figures render inline on GitHub.

---

## Quick start

```bash
git clone <repo-url>
cd scRNA_MSK_TME_Analysis
pip install -r requirements.txt
```

Download the per-sample `.h5ad` files for GSE123902 from GEO (or convert the supplied count matrices with `scanpy.read_10x_mtx`), place them in one directory, and set two paths at the top of the notebook:

```python
data_dir = r".../cache1_h5ad"   # one .h5ad per sample
out_dir  = r".../results"       # figures and DEG table are written here
```

Then run `scRNA_pipeline.ipynb` top to bottom. Runtime is a few minutes on a laptop CPU; Harmony and Leiden dominate.

---

## Limitations

These are stated up front because they bound what the figures can support.

- **Cell-level DE, not patient-level.** Wilcoxon on ~40k cells treats each cell as an independent observation, so *P*-values are inflated by orders of magnitude; the meaningful unit of replication is the 12 patients. A pseudobulk (sample-level) test is the correct primary analysis and is not included here.
- **Composition confounds expression.** Tumour and non-tumour samples differ in cell-type make-up, so a "DEG between conditions" partly reflects *which cells* were sampled rather than a state change within a cell type. Cell-type-stratified DE is needed to separate the two.
- **No malignant-cell inference.** The annotation uses an immune-focused model; epithelial cells are labelled but not split into malignant vs. normal epithelium (no `inferCNV` / `CopyKAT`). Malignant-cell-specific statements are therefore not made.
- **Heat-shock signal may be technical.** `HSPA1A`/`HSP90AA1`/`DNAJB1` induction is a well-documented artefact of tissue dissociation and hypoxia, and tumour specimens generally experience longer ischaemia. The chaperone module is reported as observed, not interpreted as biology.
- **Single dataset, no external validation.** No held-out cohort and no down-sampling control.
- **Small *n* for composition.** Four non-tumour samples, and the tumour group spans a wide range of epithelial purity; group-level composition statements are descriptive only.
- **DPT root is arbitrary** (first cell of Leiden cluster 0), so pseudotime values are not directly comparable to other studies.

---

## What I would do next

1. **Pseudobulk DE** with `decoupler` + limma-voom or DESeq2 at sample level, and report effect sizes with patient-level confidence intervals.
2. **Cell-type-stratified DE** (e.g. macrophages only, T cells only) to remove the composition confound.
3. **Malignant vs. normal epithelium** via `inferCNV` or `CopyKAT`, which would make the "tumour cell" claims in Figure 2 exact rather than approximate.
4. **Composition testing** with a proper model (scCODA / milo) instead of stacked bar eyeballing — 12 samples is enough for a Bayesian compositional test, not enough for per-sample bar charts alone.
5. **Trajectory within one lineage** (monocyte → macrophage), with the root chosen from marker expression and validated by RNA velocity.

---

## Data & citation

Data: NCBI GEO accession **GSE123902** — *The single cell transcriptional landscape of human lung adenocarcinoma (primary tumours and metastases)*. Please cite the original study when reusing the data; this repository only reanalyses it.

Code: free to reuse. If you spot an error in the analysis, an issue is very welcome.
