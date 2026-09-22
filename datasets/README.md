# Datasets

Sections that are not bundled inside the R package but are needed to reproduce the
bundled example differential-expression tables.

| file | source | used for |
|---|---|---|
| `P2N_Liver.rds` | Liu et al., J Hepatol 2023 (doi:10.1016/j.jhep.2023.01.011) — tumour-adjacent normal liver | example DEG table `02_P2N_liver_TLS_ROI_vs_rest.csv` |

Redistributed from Mendeley Data (doi:10.17632/skrx2fz79n.1) under CC BY 4.0.

This is a reduced copy of the deposited object: genes with no counts in any spot were
removed to keep the file small. Every gene that any analysis can test is retained, so
results are unchanged.

The two case-study sections ship inside the package at `inst/extdata/`
(`case_study1_CRLM.rds`, `case_study2_OSCC_spaceranger.zip`), and the ROI spot indices
behind every example table are at `inst/extdata/example_rois/` and in
[`roi_indices/`](roi_indices/).

To reproduce table 02: upload `P2N_Liver.rds`, import `roi_indices/P2N_TLS_157spots.csv`
with **Load ROI index (.csv)** on the map, then run that ROI versus Rest.
