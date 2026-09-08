# Datasets

Sections that are not bundled inside the R package but are needed to reproduce the
bundled example differential-expression tables.

| file | source | used for |
|---|---|---|
| `P2N_Spatial.rds` | Liu et al., J Hepatol 2023 (doi:10.1016/j.jhep.2023.01.011) — tumour-adjacent normal liver | example DEG table `02_P2N_liver_TLS_ROI_vs_rest.csv` |

The two case-study sections ship inside the package at `inst/extdata/`
(`case_study1_CRLM.rds`, `case_study2_OSCC.rds`), and the ROI spot indices behind every
example table are at `inst/extdata/example_rois/`.

To reproduce table 02: upload `P2N_Spatial.rds`, choose its ROI under
"Load ROI (.csv)" on the map, then run that ROI versus Rest.
