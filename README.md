# SpatialROI <img src="man/figures/logo_spatialroi.png" align="right" height="130" />

SpatialROI is designed to facilitate the ROI-specific exploration, visualization, and analysis of spatial transcriptomics data.

**SpatialROI** is an interactive R package with a browser-based interface that enables spatial visualization and analysis directly from processed Seurat objects or SpaceRanger output. Users can interactively select regions of interest (ROI), visualize gene expression patterns, and perform downstream analyses such as cell-type signature scoring, clustering, and differential expression analysis. In addition to the GUI, SpatialROI also provides a set of modular functions for scripted workflows, enabling customized analyses. Usage examples for these functions are provided in the [Function Workflow Vignette](vignettes_functions.Rmd).

SpatialROI can be accessed via a public demo hosted by the University of Pittsburgh: [https://shiny.crc.pitt.edu/spatial_api/](https://shiny.crc.pitt.edu/spatial_api/).

**Hosted-use note:** The public instance is intended for demonstration. For
unpublished, patient-derived, large, or multi-section data, use a local installation.

---

## Installation

```r
# Bioconductor packages (Seurat, GSVA)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Seurat", "GSVA"))

# Installs SpatialROI along with optional dependencies (spacexr, presto)
devtools::install_github("myaol/SpatialROI", dependencies = TRUE)
```

---

## Launching the App

### With Example Data

```r
library(SpatialROI)
run_spatial_selector("demo")
```

### With Your Own Data

```r
library(SpatialROI)
Seurat_object <- readRDS("path/to/your_seurat.rds")
# Update the object to current Seurat version
Seurat_object <- UpdateSeuratObject(Seurat_object)
run_spatial_selector(Seurat_object, sample_name = "MyExperiment", show_image = TRUE)
```

Or load data from the **Visualization** panel of the app, under **Data Source**:

- **10x Visium Seurat object (.rds)**: a processed object with a spatial image, spot
  coordinates and log-normalised expression.
- **10x Visium Space Ranger output (.zip)**: the Space Ranger output folder (filtered
  feature-barcode `.h5` plus the `spatial/` folder, with no `.gz` files), compressed as one
  `.zip`. Spots are quality-filtered and log-normalised on load.

### Supported input and size

- SpatialROI is designed and tested for **10x Genomics Visium**. A Seurat object
  must retain a Visium spatial image and tissue coordinates. Raw input must be a
  SpaceRanger output bundle containing the filtered feature-barcode matrix and
  `spatial/` files.
- **Other spatial transcriptomics platforms are not currently supported.** Seurat
  can represent data from many spatial technologies, but representation in a
  Seurat object does not imply compatibility with SpatialROI: imaging-based
  platforms (Xenium, CosMx, MERSCOPE), Slide-seq, and Visium HD bin structures
  have different data structures and are not validated here. The app stops with
  an explicit message when the uploaded spatial image is not recognized as
  Visium.
- Uploads are limited to **150 MB on the public server** and **500 MB in a local
  installation**. Use a local installation for large objects.
- RDS files expand in memory. The 32-MB example requires approximately 430 MB after
  loading; memory requirements increase with spots, assays, and image size.

### Bundled example dataset

The example is one human colorectal cancer 10x Visium tissue section from
Valdeolivas et al. (2024), with **17,529 genes and 1,253 tissue spots**, an H&E
image, SCT-normalized expression, and precomputed broad-cell-type proportions.
The associated publication is [*npj Precision Oncology* 8, 7
(2024)](https://doi.org/10.1038/s41698-023-00488-4).

---

## ROI Selection Tool

### Quick Spot Selection with `draw_ROI()`

If you only need to select spots from a region of interest without launching the full analysis app:

```r
library(SpatialROI)

# Load your Seurat object
Seurat_object <- readRDS("path/to/your_seurat.rds")

# Launch interactive ROI selector
selected_spots <- draw_ROI(Seurat_object, sample_name = "MyExperiment")

# The function returns a vector of spot IDs
print(selected_spots)
length(selected_spots)

# Use the selected spots for downstream analysis
subset_data <- subset(Seurat_object, cells = selected_spots)
```

This function supports multiple ROI selections and returns a vector of spot IDs, ideal for custom downstream workflows.

---

## Features

- 🗺️ **ROI Drawing** - Freehand drawing tools to select custom regions of interest
- 🧬 **Gene Set and Pathway Visualization** - Spatially map custom gene lists, cell-type signatures, or pathway gene sets
- 🔗 **Ligand-Receptor Colocalization** - ROI-specific, Gaussian-smoothed ligand-receptor score analysis
- 🧩 **Cell Type Deconvolution** - RCTD-based cell type deconvolution within user-defined ROIs
- 📊 **Spot Clustering** - Identify spatial domains using graph-based clustering methods
- 📈 **DEG Analysis** - Find differentially expressed genes between selected groups or clusters
- ⚖️ **Feature Comparison** - Statistical comparison plots with parametric/non-parametric tests
- 🧩 **Multi-Sample Comparison** - Compare ROI-versus-rest DEG tables across sections or studies without pooling expression data
- 💾 **Data Export** - Save and reload ROI spot indices, and download DEG results and Seurat subsets
- 🖼️ **Figure Export** - Download UMAP, volcano, Moran, violin, and heatmap figures as PDFs

---

## Example Data and ROI Spot Indices

Three datasets ship inside the package and load from the interface without any
upload:

| button | dataset | source |
|---|---|---|
| Default Data (CRC) | colorectal cancer Visium, 1,253 spots | Valdeolivas *et al.* 2024 |
| Case Study 1 (CRLM) | colorectal cancer liver metastasis, 3,721 spots | Wu *et al.* 2022, OEP00001756 |
| Case Study 2 (OSCC) | oral squamous cell carcinoma, raw Space Ranger output | Arora *et al.* 2023, GSE208253 |

### Reproducing the Multi-Sample example tables

The Multi-Sample panel ships three ROI-versus-rest differential-expression
tables. The exact spot barcodes behind each one are published here:

**[`datasets/roi_indices/`](datasets/roi_indices/)**

| ROI index | spots | dataset to load | regenerates |
|---|---|---|---|
| `CRC_TLS_41spots.csv` | 41 | Default Data (CRC) — bundled | `01_CRC_TLS_ROI_vs_rest.csv` |
| `CRLM_TLS_86spots.csv` | 86 | Case Study 1 (CRLM) — bundled | `03_CRLM_liver_TLS_ROI_vs_rest.csv` |
| `P2N_TLS_157spots.csv` | 157 | [`datasets/P2N_Liver.rds`](datasets/) | `02_P2N_liver_TLS_ROI_vs_rest.csv` |

To reproduce a table: load the dataset, import its index with **Load ROI index
(.csv)** on the map, then run that ROI versus Rest. `CRC_TLS_41spots.csv` is also
the 41-spot region used in the cross-tool output-concordance analysis.

Any `.csv` with a `spot_id` column imports, so regions defined in other software
can be loaded the same way; the region takes the file's name unless the file
carries a `roi` column.

### Upload limits

Uploads are limited to 150 MB on the public server and 500 MB in a local
installation. A Seurat object needs roughly three times its file size in memory
once loaded, so large sections are better analysed locally. The local limit can
be raised before launching:

```r
options(shiny.maxRequestSize = 2 * 1024^3)   # 2 GB
```

## Reference Datasets

Two RCTD references are built into SpatialROI and can be selected in the deconvolution
panel: colorectal cancer (CRC, from GSE132465) and colorectal cancer liver metastasis
(CRLM, from GSE225857). Curated references for lung adenocarcinoma, lung squamous cell
carcinoma, renal cell carcinoma, breast cancer, hepatocellular carcinoma, oral squamous
cell carcinoma, CRLM and mouse brain are hosted on Zenodo:

DOI: https://doi.org/10.5281/zenodo.20554051

These datasets can be downloaded separately and supplied to SpatialROI for RCTD-based cell-type deconvolution.
Uploaded Seurat references must contain original RNA counts and cell-type labels
in active identities or a metadata column; retained cell types require at least
25 cells. Spatial objects used to rerun RCTD must contain an original `Spatial`
or `RNA` raw-count assay.

---

## Documentation

📚 **Detailed tutorials and examples:**
- [User Guide Vignette](vignettes.Rmd) - GUI Step-by-step walkthrough
- [Function Workflow Vignette](vignettes_functions.Rmd) - Scripted workflow
- [Manuscript](https://academic.oup.com/bioinformaticsadvances) - Lu et al. 2026, *Bioinformatics Advances* (in submission)

Static statistical plots are exported as publication-ready PDFs.



---

## Limitations and Future Work

- **Platform scope.** SpatialROI currently supports 10x Genomics Visium.
  Extending support to additional spatial transcriptomics platforms, including
  Visium HD, Xenium, CosMx, MERSCOPE, and Slide-seq, is a potential direction
  for future development and would require platform-specific implementation and
  validation.
- **Single-section analysis.** Interactive ROI drawing and most downstream
  analyses are performed on one image-aligned tissue section at a time. The
  Multi-Sample module supports cross-section comparison using ROI-versus-rest
  DEG results without pooling raw expression matrices.
- **Cross-sample interpretation.** Multi-Sample comparisons are descriptive and
  may be influenced by batch effects, patient-level differences, tissue
  composition, and ROI annotation.
- **Statistical interpretation.** Spot-level statistical tests describe
  within-section differences and should not be interpreted as independent
  patient-level replication.
- **Large datasets.** Local use is recommended for large or unpublished
  datasets, particularly when server upload or memory limits may become
  restrictive.

---

## Getting Help

If you encounter bugs or have suggestions for improvements:
- **Report issues:** [GitHub Issues](https://github.com/myaol/SpatialROI/issues)
- **Contact authors:**
  - Mengyao Lu: [mel373@pitt.edu](mailto:mel373@pitt.edu)
  - Aodong Qiu: [qiuaodon@pitt.edu](mailto:qiuaodon@pitt.edu)
  - Lujia Chen: [luc17@pitt.edu](mailto:luc17@pitt.edu)

When reporting issues, please include your sessionInfo(), a minimal reproducible example, and any error messages.

---

## Citation

If you use SpatialROI in your research, please cite:

```bibtex
@article{SpatialROI2026,
  title = {SpatialROI: An Interactive R Shiny Package for Manual Region-Based Analysis of Spatial Transcriptomics Data},
  author = {Lu, Mengyao and Qiu, Aodong and Lu, Xinghua and Xu, Min and Chen, Lujia},
  journal = {Bioinformatics Advances},
  year = {2026},
  note = {manuscript in submission},
  url = {https://github.com/myaol/SpatialROI}
}
```

---

## Disclaimer

SpatialROI is designed to facilitate intuitive visualization, region selection, and exploratory analysis of spatial transcriptomics data. The tool provides convenient interfaces for clustering, differential expression, and feature comparison, but these analyses are intended for exploratory purposes only. Users should validate any biological interpretations using appropriate statistical or experimental methods.

---

## Acknowledgments

This research was supported in part by the University of Pittsburgh Center for Research Computing and Data, RRID:SCR_022735, through the resources provided. Specifically, this work used the HTC cluster, which is supported by NIH award number S10OD028483.

This work was supported by NIH grants including NHGRI R01HG014023, NLM 4R00LM013089, 5R01LM012011, and by U.S. NIH grants R35GM158094 and R01GM134020, as well as NSF grants DBI-2238093, DBI-2422619, IIS-2211597, and MCB-2205148.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
