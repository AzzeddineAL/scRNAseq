# Glioblastoma scRNA-seq: tumour states and microenvironment

**Authors:** Azzeddine ALILICHE, Anais HLUSZKO, Corentin OSTA

Single-cell RNA-seq analysis of a glioblastoma sample (10x Genomics Chromium Flex, FFPE tissue) in R / Seurat. The question we ask: **in this tumour, which tumour states are present, what microenvironment surrounds them, and which signals might connect them?**

## Data

- 10x Genomics public dataset: *Mixture of Lung Cancer and Glioblastoma FFPE Tissues Dissociated Manually or using gentleMACS Dissociator, Multiplexed Samples, 4 Probe Barcodes* (Chromium Flex, Cell Ranger 7.1.0).
- We analyse **one sample: the glioblastoma, manual dissociation (probe barcode BC3)**. [Raw feature-barcode matrix](https://cf.10xgenomics.com/samples/cell-exp/7.1.0/4plex_human_lungcancer_glioblastoma_scFFPE_multiplex_Glioblastoma_Manual_BC3/4plex_human_lungcancer_glioblastoma_scFFPE_multiplex_Glioblastoma_Manual_BC3_count_sample_raw_feature_bc_matrix.tar.gz)

## Pipeline

| Step | Document | Methods |
|---|---|---|
| 1. Quality control | `Glioblastoma_scRNAseq_QC.qmd` | Gene / UMI / mitochondrial thresholds, doublet removal (scDblFinder) |
| 2. Clustering & annotation | `Glioblastoma_scRNAseq_Clustering.qmd` | PCA, Louvain clustering, UMAP, SingleR, GBmap label transfer, marker genes |
| 3. Malignant cells | `Glioblastoma_scRNAseq_CNV.qmd` | inferCNV, chr7 gain / chr10 loss score |
| 4. Tumour states | `Glioblastoma_scRNAseq_States.qmd` | Neftel et al. 2019 meta-modules, GO enrichment (g:Profiler) |
| 5. Microenvironment | `Glioblastoma_scRNAseq_TME.qmd` | TAM origin and sub-populations, T cell markers, GO enrichment |
| 6. Cell-cell communication | `Glioblastoma_scRNAseq_CellChat.qmd` | CellChat v2 |
| Synthesis | `Glioblastoma_scRNAseq_Synthesis.qmd` | Key results and conclusion |

The rendered reports (self-contained HTML) are in [`report/`](report/). Start with `Glioblastoma_scRNAseq_Synthesis.html`. GitHub doesn't display HTML files, so download them and open them in a browser.

## Key results

- **4005 cells** after QC, 11 cell populations: tumour cells, macrophages (TAMs), vascular cells, T/NK cells and infiltrated brain tissue (neurons, oligodendrocytes, astrocytes, OPC).
- **1268 malignant cells** confirmed by inferCNV: chr7 gain and chr10 loss, chr13 loss, and a focal chr12 gain, probably an MDM2 amplification.
- **Three main tumour states:** MES-like (40%), NPC-like (33%) and AC-like (26%). AC-like and MES-like cells form a continuum, while NPC-like cells are a distinct group that holds most of the cycling cells. **MES-like cells are hypoxic** (top GO terms: response to hypoxia, angiogenesis).
- **TAMs:** mostly bone marrow-derived (63%) and immunosuppressive, with microglia-like, hypoxic and MRC1+ sub-populations. **T cells are rare and poorly active**, a "cold" immune environment.
- **Communication:** SPP1 → CD44 (TAMs → MES-like / AC-like cells) and VEGFA → FLT1 (MES-like cells and hypoxic TAMs → vascular cells) are supported. OSM → OSMR and PD-L1 → PD-1 are not detected.

**Model for this tumour:** a hypoxic niche in which MES-like tumour cells and hypoxic, immunosuppressive TAMs interact through SPP1–CD44 and both drive angiogenesis through VEGFA.

**Limits:** a single tumour with no replicate, so the results are descriptive and generate hypotheses. All results are inferred from RNA, and the Flex chemistry only measures genes in its probe panel.

## Repository structure

```
code/                 Quarto documents (one per step), _quarto.yml, prep_gbmap_ref.R
data/raw_data/        10x raw feature-barcode matrix
data/reference/       inferCNV gene positions (hg38), Neftel et al. Table S2
data/processed_data/  Seurat objects saved by each step
report/               Rendered HTML reports
```

## Reproduce the analysis

Requirements: R 4.5, Quarto (bundled with RStudio) and the R packages Seurat, scDblFinder, SingleR, celldex, infercnv, CellChat, gprofiler2, readxl, rhdf5.

On Windows you also need:

- **Rtools 4.5**, to build CellChat from GitHub (`pak::pak("jinworks/CellChat")`).
- **JAGS 4.x**, for rjags and inferCNV. rjags doesn't work with JAGS 5 yet.
- **Python 3**, while installing inferCNV only, through the `PYTHON` environment variable if R doesn't find it.

Steps:

1. Download the **GBmap core atlas** (`.h5ad`, 8 GB) from [CELLxGENE](https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c) into `data/gbmap_core/`. Then build the light reference once, from the `code/` folder:
   ```
   Rscript prep_gbmap_ref.R
   ```
2. Render the whole project in pipeline order, from the `code/` folder:
   ```
   quarto render
   ```
   Or use *Build → Render Project* in RStudio. The full render takes about 20 minutes.

GO enrichment calls the g:Profiler web service, so an internet connection is needed.

## References

- Neftel C. *et al.* An integrative model of cellular states, plasticity, and genetics for glioblastoma. *Cell* 178, 835–849 (2019).
- Hara T. *et al.* Interactions between cancer cells and immune cells drive transitions to mesenchymal-like states in glioblastoma. *Cancer Cell* 39, 779–792 (2021).
- Ruiz-Moreno C. *et al.* Charting the single-cell and spatial landscape of IDH-wild-type glioblastoma with GBmap. *Neuro-Oncology* (2025).
- Jin S. *et al.* CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics. *Nature Protocols* (2024).
- Hao Y. *et al.* Dictionary learning for integrative, multimodal and scalable single-cell analysis (Seurat v5). *Nature Biotechnology* (2024).
- inferCNV of the Trinity CTAT Project, Broad Institute.
