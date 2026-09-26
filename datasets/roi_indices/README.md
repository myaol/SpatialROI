# ROI spot-index files

The exact spot barcodes behind each bundled example DEG table. Download one, load the
matching dataset in SpatialROI, then import it with **⬆ Load ROI index (.csv)** on the map
and run that ROI versus Rest to regenerate the table.

| file | spots | dataset to load | regenerates |
|---|---|---|---|
| `CRC_TLS_41spots.csv` | 41 | Default Data (CRC) — bundled | `01_CRC_TLS_ROI_vs_rest.csv` |
| `CRLM_TLS_86spots.csv` | 86 | Case Study 1 (CRLM) — bundled | `03_CRLM_liver_TLS_ROI_vs_rest.csv` |
| `P2N_TLS_157spots.csv` | 157 | `../P2N_Liver.rds` — download first | `02_P2N_liver_TLS_ROI_vs_rest.csv` |

`CRC_TLS_41spots.csv` is also the 41-spot region used in the cross-tool output-concordance
benchmark.

## Case Study 1

`case_study1/` holds the four regions used in Case Study 1 of the supplementary material.
Load **Case Study 1 (CRLM)**, then import each file with **⬆ Load ROI index (.csv)**. The
regions keep their names (ROI 1 to ROI 4).

| file | spots |
|---|---|
| `case_study1/CaseStudy1_CRLM_ROI_1_spot_index.csv` | 90 |
| `case_study1/CaseStudy1_CRLM_ROI_2_spot_index.csv` | 62 |
| `case_study1/CaseStudy1_CRLM_ROI_3_spot_index.csv` | 65 |
| `case_study1/CaseStudy1_CRLM_ROI_4_spot_index.csv` | 59 |

Format: `sample`, `roi`, `spot_id`. Any file with a `spot_id` column will import, so
regions defined in other software can be loaded the same way.
