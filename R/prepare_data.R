#' Prepare Seurat data for spatial selector
#' @param seurat_input Seurat object or data list
#' @param sample_name Character label for dataset
#' @param show_image Logical; whether to show image
#' @return A list containing processed Seurat object and metadata
#' @keywords internal
#' @importFrom sf st_as_sf st_geometry
#' @importFrom magrittr %>%


prepare_seurat_data <- function(seurat_input, sample_name = "sample", show_image = TRUE) {

  # 1. Normalize input
  if (is.list(seurat_input) && "seurat" %in% names(seurat_input)) {
    seurat_obj <- seurat_input$seurat
  } else if (inherits(seurat_input, "Seurat")) {
    seurat_obj <- seurat_input
  } else {
    stop("Input must be either a Seurat object or a list containing a $seurat component.")
  }

  # 2. Sanity checks
  if (length(seurat_obj@images) == 0) {
    stop("No spatial images found in the Seurat object.")
  }

  # 3. Extract image + coordinates safely
  image_name <- names(seurat_obj@images)[1]
  coords <- NULL

  # Try the modern Seurat helper
  suppressWarnings({
    try({
      coords <- Seurat::GetTissueCoordinates(seurat_obj, image = image_name)
    }, silent = TRUE)
  })

  # Try legacy slot
  if (is.null(coords)) {
    try({
      coords <- seurat_obj@images[[image_name]]@coordinates
    }, silent = TRUE)
  }

  # Try if coordinates are stored in metadata
  if (is.null(coords)) {
    try({
      if (all(c("imagerow", "imagecol") %in% colnames(seurat_obj@meta.data))) {
        coords <- seurat_obj@meta.data[, c("imagerow", "imagecol"), drop = FALSE]
      }
    }, silent = TRUE)
  }

  # Stop if nothing worked
  if (is.null(coords) || nrow(coords) == 0) {
    stop("No spatial coordinates found in seurat_input — please check that it contains spatial data.")
  }

  # Process coordinates
  spots_df <- data.frame(spot_id = rownames(coords), stringsAsFactors = FALSE)
  if ("imagerow" %in% colnames(coords) && "imagecol" %in% colnames(coords)) {
    spots_df$x <- coords$imagecol
    spots_df$y <- coords$imagerow
  } else if ("row" %in% colnames(coords) && "col" %in% colnames(coords)) {
    spots_df$x <- coords$col
    spots_df$y <- coords$row
  } else {
    spots_df$x <- coords[,1]
    spots_df$y <- coords[,2]
  }

  spots_sf <- sf::st_as_sf(spots_df, coords = c("x","y"), crs = NA)
  coords_matrix <- as.matrix(do.call(rbind, sf::st_geometry(spots_sf)))

  spots_sf$x <- coords_matrix[, 1]
  spots_sf$y <- coords_matrix[, 2]
  spots_sf$y <- max(spots_sf$y) - spots_sf$y + min(spots_sf$y)

  # H&E image processing
  he_image_base64 <- NULL
  he_image_bounds <- NULL

  if (show_image) {
    tryCatch({
      image_obj <- seurat_obj@images[[image_name]]
      if (class(image_obj)[1] == "VisiumV1") {
        he_image_data <- image_obj@image
        coords_full <- Seurat::GetTissueCoordinates(seurat_obj, image = image_name)

        if ("pxl_col_in_fullres" %in% colnames(coords_full)) {
          pixel_x <- coords_full$pxl_col_in_fullres
          pixel_y <- coords_full$pxl_row_in_fullres
        } else {
          pixel_x <- spots_sf$x
          pixel_y <- spots_sf$y
        }

        pixel_x_min <- min(pixel_x, na.rm = TRUE)
        pixel_x_max <- max(pixel_x, na.rm = TRUE)
        pixel_y_min <- min(pixel_y, na.rm = TRUE)
        pixel_y_max <- max(pixel_y, na.rm = TRUE)
        x_buffer_px <- (pixel_x_max - pixel_x_min) * 0.1
        y_buffer_px <- (pixel_y_max - pixel_y_min) * 0.1

        crop_x_min <- max(1, floor(pixel_x_min - x_buffer_px))
        crop_x_max <- min(dim(he_image_data)[2], ceiling(pixel_x_max + x_buffer_px))
        crop_y_min <- max(1, floor(pixel_y_min - y_buffer_px))
        crop_y_max <- min(dim(he_image_data)[1], ceiling(pixel_y_max + y_buffer_px))

        he_image_cropped <- he_image_data[crop_y_min:crop_y_max, crop_x_min:crop_x_max, ]

        temp_file <- tempfile(fileext = ".png")
        png(temp_file, width = dim(he_image_cropped)[2], height = dim(he_image_cropped)[1])
        par(mar = c(0,0,0,0))
        plot(as.raster(he_image_cropped))
        dev.off()

        he_image_base64 <- paste0("data:image/png;base64,", base64enc::base64encode(temp_file))
        unlink(temp_file)

        he_image_bounds <- list(
          south = crop_y_max,
          west = crop_x_min,
          north = crop_y_min,
          east = crop_x_max
        )
      }
    }, error = function(e) {
      message("Could not extract H&E image: ", e$message)
    })
  }

  # =====================
  # Enhanced Signature Library - HUMAN
  # Cell marker gene signatures curated from CellMarker 2.0 database and literature
  # Citation: Hu C, Li T, Xu Y, et al. CellMarker 2.0: an updated database of
  # manually curated cell markers in human/mouse and web tools based on scRNA-seq data.
  # Nucleic Acids Res. 2023;51(D1):D870-D876. doi: 10.1093/nar/gkac947
  # Database: http://bio-bigdata.hrbmu.edu.cn/CellMarker/
  # Download: http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html
  #
  # CellMarker 2.0 contains 26,915 cell markers, 2,578 cell types, and 656 tissues
  # For comprehensive cell marker information, please visit the database website

  signature_library_human <- list(
    "Custom" = list(genes = c()),

    # ========== IMMUNE CELL SIGNATURES - EXPANDED ==========

    # B cell lineage
    "Immune: B cells" = list(
      genes = c("CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "CD22", "BANK1",
                "BLK", "FCRL5", "VPREB3", "IGHM", "IGHD")
    ),

    "Immune: B cells naive" = list(
      genes = c("IGHD", "IGHM", "CCR7", "SELL", "CD19", "MS4A1", "TCL1A",
                "FCER2", "IL4R")
    ),

    "Immune: B cells memory" = list(
      genes = c("CD27", "CD19", "MS4A1", "TNFRSF13B", "CD80", "CD73", "AIM2")
    ),

    "Immune: B cells activated" = list(
      genes = c("CD69", "CD83", "JUN", "IRF4", "MYC", "CD86", "TNFRSF17")
    ),

    "Immune: Plasma cells" = list(
      genes = c("MZB1", "SDC1", "CD38", "XBP1", "PRDM1", "JCHAIN", "IGHG1",
                "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGKC", "IGLC2",
                "IGLC3", "SSR4", "DERL3", "FKBP11", "PRDX4", "SLAMF7", "CD27")
    ),

    "Immune: Plasmablasts" = list(
      genes = c("IGHG1", "IGHG4", "JCHAIN", "IGLC3", "MZB1", "IGKC", "IGHM",
                "XBP1", "CD38", "CD27", "SDC1")
    ),

    # T cell lineage
    "Immune: T cells CD8" = list(
      genes = c("CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "GZMK", "CCL5",
                "GZMA", "PRF1", "NKG7", "TRAC", "TRBC1", "TRBC2")
    ),

    "Immune: T cells CD8 naive" = list(
      genes = c("CD8A", "CD8B", "CCR7", "SELL", "LEF1", "TCF7", "CD27",
                "CD28", "IL7R")
    ),

    "Immune: T cells CD8 effector" = list(
      genes = c("CD8A", "GZMB", "GZMH", "PRF1", "GNLY", "NKG7", "FGFBP2",
                "FCGR3A", "CX3CR1")
    ),

    "Immune: T cells CD8 exhausted" = list(
      genes = c("CD8A", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX",
                "ENTPD1", "LAYN", "CXCL13")
    ),

    "Immune: T cells CD4" = list(
      genes = c("CD4", "IL7R", "CD3D", "CD3E", "CD3G", "CD40LG", "ICOS",
                "MAF", "TRAC", "TRBC1", "TRBC2", "CD28")
    ),

    "Immune: T cells CD4 naive" = list(
      genes = c("CD4", "CCR7", "SELL", "LEF1", "TCF7", "CD27", "CD28",
                "IL7R", "MAL")
    ),

    "Immune: T cells CD4 memory" = list(
      genes = c("CD4", "IL7R", "CCR7", "SELL", "CD27", "CD28", "TNFRSF4",
                "LTB", "AQP3")
    ),

    "Immune: T cells regulatory (Tregs)" = list(
      genes = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TNFRSF4",
                "IL10", "TIGIT", "LAYN", "CCR8", "IL2RB", "CD4")
    ),

    "Immune: T cells follicular helper (Tfh)" = list(
      genes = c("CXCR5", "PDCD1", "BCL6", "CD4", "IL21", "ICOS", "CD40LG",
                "MAF", "TOX2", "BTLA")
    ),

    "Immune: T cells Th1" = list(
      genes = c("CD4", "IFNG", "TBX21", "IL12RB2", "CXCR3", "CCR5", "TNF")
    ),

    "Immune: T cells Th2" = list(
      genes = c("CD4", "GATA3", "IL4", "IL5", "IL13", "CCR4", "IL1RL1")
    ),

    "Immune: T cells Th17" = list(
      genes = c("CD4", "RORC", "IL17A", "IL17F", "IL23R", "CCR6", "IL22")
    ),

    # NK cells
    "Immune: NK cells" = list(
      genes = c("NCAM1", "FCGR3A", "NKG7", "GNLY", "KLRD1", "KLRF1", "KLRB1",
                "KLRC1", "KLRC2", "KLRC3", "NCR1", "GZMB", "PRF1", "XCL1",
                "XCL2", "FCER1G", "TYROBP")
    ),

    "Immune: NK cells CD56bright" = list(
      genes = c("NCAM1", "SELL", "IL7R", "CD44", "XCL1", "XCL2", "GZMK")
    ),

    "Immune: NK cells CD56dim" = list(
      genes = c("FCGR3A", "KLRF1", "KLRD1", "PRF1", "GZMB", "S1PR5", "CX3CR1")
    ),

    # Myeloid lineage
    "Immune: Monocytes" = list(
      genes = c("CD14", "FCGR3A", "S100A8", "S100A9", "LYZ", "VCAN", "FCN1",
                "S100A12", "CSF1R", "CD68", "ITGAM", "ITGAX")
    ),

    "Immune: Monocytes classical" = list(
      genes = c("CD14", "S100A8", "S100A9", "S100A12", "LYZ", "FCN1", "VCAN",
                "ITGAM")
    ),

    "Immune: Monocytes non-classical" = list(
      genes = c("FCGR3A", "CD14", "CDKN1C", "MS4A7", "ASAH1", "IFITM3", "LST1")
    ),

    "Immune: Monocytes intermediate" = list(
      genes = c("CD14", "FCGR3A", "HLA-DRA", "HLA-DRB1", "CD74", "CCR2")
    ),

    "Immune: Macrophages" = list(
      genes = c("CD68", "CD163", "CSF1R", "MSR1", "C1QA", "C1QB", "C1QC",
                "MARCO", "MRC1", "CD14", "FCGR1A", "FCGR2A")
    ),

    "Immune: Macrophages M1" = list(
      genes = c("CD68", "CD80", "CD86", "IL1B", "IL6", "TNF", "CXCL9",
                "CXCL10", "CXCL11", "NOS2", "IDO1", "SOCS3")
    ),

    "Immune: Macrophages M2" = list(
      genes = c("CD68", "CD163", "MRC1", "MSR1", "LYVE1", "MAF", "FOLR2",
                "F13A1", "SEPP1", "STAB1", "CD209", "CLEC10A")
    ),

    "Immune: Macrophages alveolar" = list(
      genes = c("MARCO", "FABP4", "PPARG", "MCEMP1", "CD68", "SLC7A8",
                "CHIT1", "MSR1")
    ),

    "Immune: Macrophages tumor-associated (TAM)" = list(
      genes = c("CD68", "CD163", "MRC1", "MSR1", "TREM2", "SPP1", "APOE",
                "C1QA", "C1QB", "FOLR2", "LYVE1")
    ),

    "Immune: Dendritic cells" = list(
      genes = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD1C",
                "FCER1A", "CLEC10A", "CD209", "ITGAX", "ITGAM")
    ),

    "Immune: Dendritic cells conventional (cDC1)" = list(
      genes = c("CLEC9A", "XCR1", "CADM1", "BATF3", "IRF8", "IDO1", "ITGAE")
    ),

    "Immune: Dendritic cells conventional (cDC2)" = list(
      genes = c("CD1C", "FCER1A", "CLEC10A", "CD209", "SIRPA", "IRF4")
    ),

    "Immune: Dendritic cells plasmacytoid (pDC)" = list(
      genes = c("IL3RA", "CLEC4C", "NRP1", "TCF4", "IRF7", "IRF8", "LILRA4",
                "GZMB", "ITM2C")
    ),

    "Immune: Dendritic cells activated" = list(
      genes = c("CD83", "LAMP3", "CCR7", "FSCN1", "CCL19", "CCL22", "IL1B",
                "MARCKSL1", "IDO1")
    ),

    "Immune: Neutrophils" = list(
      genes = c("FCGR3B", "CSF3R", "S100A8", "S100A9", "S100A12", "CXCR2",
                "FPR1", "FCGR2A", "MPO", "ELANE", "CTSG", "PRTN3", "AZU1",
                "LCN2", "MMP8", "MMP9", "CAMP", "LTF", "DEFA3", "DEFA4")
    ),

    "Immune: Eosinophils" = list(
      genes = c("SIGLEC8", "IL5RA", "CCR3", "CLC", "EPX", "PRG2", "PRG3",
                "RNASE2", "RNASE3", "ALOX15", "IL1RL1")
    ),

    "Immune: Basophils" = list(
      genes = c("ENPP3", "HDC", "MS4A2", "GATA2", "CPA3", "TPSAB1", "IL3RA",
                "SLC18A2", "IL4", "IL13")
    ),

    "Immune: Mast cells" = list(
      genes = c("TPSAB1", "TPSB2", "CPA3", "KIT", "GATA2", "HDC", "MS4A2",
                "HPGDS", "CTSG", "RGS13", "IL1RL1", "FCER1A", "SLC18A2")
    ),

    # ========== FUNCTIONAL IMMUNE SIGNATURES - EXPANDED ==========

    "Immune: Cytotoxicity" = list(
      genes = c("PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "GNLY", "NKG7",
                "KLRK1", "KLRB1", "FASLG", "IFNG", "TNF")
    ),

    "Immune: Exhaustion" = list(
      genes = c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "TOX2",
                "ENTPD1", "ITGAE", "LAYN", "CXCL13", "CD244", "CD160",
                "BTLA", "VSIR", "TNFRSF9")
    ),

    "Immune: Inflammation" = list(
      genes = c("IL1B", "IL6", "TNF", "CXCL8", "IL1A", "PTGS2", "CCL2",
                "CCL3", "CCL4", "CXCL1", "CXCL2", "CXCL3", "CSF2", "CSF3")
    ),

    "Immune: Interferon response" = list(
      genes = c("IFITM1", "IFITM2", "IFITM3", "ISG15", "ISG20", "IFI6",
                "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
                "MX1", "MX2", "OAS1", "OAS2", "OAS3", "IRF7", "STAT1", "STAT2")
    ),

    "Immune: Complement activation" = list(
      genes = c("C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C4A",
                "C4B", "C5", "CFB", "CFD", "CFH", "CFI")
    ),

    "Immune: Phagocytosis" = list(
      genes = c("FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B",
                "MSR1", "CD36", "MARCO", "SCARB1", "MERTK", "AXL", "TYRO3")
    ),

    "Immune: Antigen presentation" = list(
      genes = c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DPA1",
                "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "HLA-DMB",
                "CD74", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9")
    ),

    "Immune: Chemotaxis" = list(
      genes = c("CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7",
                "CCR9", "CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5", "CXCR6",
                "CX3CR1", "ACKR1", "ACKR2")
    ),

    # ========== CANCER SIGNATURES - EXPANDED ==========

    "Cancer: Epithelial" = list(
      genes = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "CLDN3", "CLDN4",
                "CLDN7", "TACSTD2", "MUC1", "ELF3", "KRT7", "KRT20", "GRHL2")
    ),

    "Cancer: Proliferation" = list(
      genes = c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CCNB2", "CCNE1",
                "CCNE2", "CDC20", "AURKA", "AURKB", "PLK1", "BUB1", "BUB1B",
                "CENPE", "CENPF", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
                "MCM7", "TYMS", "TK1", "RRM2", "PTTG1", "UBE2C", "BIRC5")
    ),

    "Cancer: Hypoxia" = list(
      genes = c("HIF1A", "VEGFA", "CA9", "SLC2A1", "LDHA", "PGK1", "BNIP3",
                "BNIP3L", "PDK1", "ENO1", "PFKL", "PFKP", "SLC16A3", "NDRG1",
                "EGLN3", "LOX", "P4HA1", "P4HA2", "ADM", "ANGPTL4", "VHL")
    ),

    "Cancer: EMT (Epithelial-Mesenchymal Transition)" = list(
      genes = c("VIM", "CDH2", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1",
                "ZEB2", "FN1", "MMP2", "MMP3", "MMP9", "MMP14", "TCF4",
                "FOXC2", "SPARC", "CDH11", "THY1", "ITGA5", "ITGB1")
    ),

    "Cancer: Stemness" = list(
      genes = c("CD44", "ALDH1A1", "PROM1", "SOX2", "NANOG", "POU5F1",
                "KLF4", "MYC", "LGR5", "BMI1", "NOTCH1", "NOTCH2", "JAG1",
                "DLL1", "CD24", "EPCAM", "ITGA6", "THY1")
    ),

    "Cancer: Angiogenesis" = list(
      genes = c("VEGFA", "VEGFB", "VEGFC", "KDR", "FLT1", "FLT4", "PECAM1",
                "ANGPT1", "ANGPT2", "TEK", "NRP1", "NRP2", "PGF", "HIF1A",
                "PDGFA", "PDGFB", "PDGFRA", "PDGFRB", "FGF2", "FGFR1")
    ),

    "Cancer: Metastasis" = list(
      genes = c("MMP2", "MMP9", "MMP11", "MMP14", "CXCR4", "CXCL12", "TWIST1",
                "SNAI1", "ZEB1", "VIM", "FN1", "CDH2", "ITGA5", "ITGB3",
                "ITGAV", "CD44", "MET", "PTGS2", "CTNNB1")
    ),

    "Cancer: Apoptosis resistance" = list(
      genes = c("BCL2", "BCL2L1", "MCL1", "BIRC5", "XIAP", "BIRC2", "BIRC3",
                "MDM2", "TP53", "BAX", "BAK1", "BID", "CASP3", "CASP8",
                "CASP9", "CYCS")
    ),

    "Cancer: DNA damage response" = list(
      genes = c("BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "RAD51",
                "XRCC1", "PARP1", "H2AX", "TP53BP1", "MDC1", "FANCD2",
                "FANCA", "FANCG")
    ),

    # ========== STROMAL SIGNATURES - EXPANDED ==========

    "Stroma: Fibroblasts" = list(
      genes = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "THY1",
                "VIM", "ACTA2", "TAGLN", "MYL9", "FN1", "SPARC")
    ),

    "Stroma: Fibroblasts (CAF)" = list(
      genes = c("FAP", "PDGFRA", "PDPN", "COL1A1", "COL1A2", "COL3A1",
                "COL5A1", "COL6A1", "FN1", "DCN", "LUM", "ACTA2", "TAGLN",
                "TGFB1", "TGFB2", "CXCL12", "CXCL14", "IL6", "CCL2",
                "POSTN", "THBS2", "MMP11", "LRRC15", "COMP")
    ),

    "Stroma: Myofibroblasts" = list(
      genes = c("ACTA2", "TAGLN", "MYL9", "MYH11", "CNN1", "ACTG2", "DES",
                "PDGFRB", "RGS5")
    ),

    "Stroma: Endothelial" = list(
      genes = c("PECAM1", "VWF", "CDH5", "KDR", "FLT1", "ENG", "MCAM",
                "PLVAP", "TEK", "ESAM", "CLDN5", "CD34", "ERG", "FLI1",
                "ICAM2", "EMCN", "NOS3", "RAMP2", "RAMP3")
    ),

    "Stroma: Endothelial arterial" = list(
      genes = c("GJA5", "BMX", "EFNB2", "SOX17", "DLL4", "HEY1", "SEMA3G")
    ),

    "Stroma: Endothelial venous" = list(
      genes = c("NR2F2", "ACKR1", "SELP", "VWF", "CLU", "APLNR")
    ),

    "Stroma: Endothelial lymphatic" = list(
      genes = c("PROX1", "LYVE1", "PDPN", "FLT4", "CCL21", "TFF3", "MMRN1")
    ),

    "Stroma: Endothelial tip cells" = list(
      genes = c("APLN", "ESM1", "ANGPT2", "FLT1", "APLN", "PGF", "KDR",
                "PDGFB", "DLL4")
    ),

    "Stroma: Pericytes" = list(
      genes = c("PDGFRB", "RGS5", "ACTA2", "MCAM", "NOTCH3", "CSPG4",
                "ABCC9", "KCNJ8", "NDUFA4L2", "COX4I2", "STEAP4", "ANGPT1")
    ),

    "Stroma: Smooth muscle cells" = list(
      genes = c("ACTA2", "MYH11", "TAGLN", "CNN1", "ACTG2", "MYL9", "DES",
                "MYLK", "SMTN", "TPM2", "TPM1", "LMOD1")
    ),

    "Stroma: Adipocytes" = list(
      genes = c("ADIPOQ", "LEP", "PLIN1", "FABP4", "LPL", "PPARG", "CEBPA",
                "CIDEC", "PLIN4", "LIPE", "DGAT2", "SLC2A4")
    ),

    # ========== TISSUE-SPECIFIC SIGNATURES ==========

    "Tissue: Neural" = list(
      genes = c("TUBB3", "MAP2", "SYP", "SNAP25", "SYT1", "NEFL", "ENO2",
                "DCX", "NCAM1", "GAP43", "STMN2", "RBFOX3")
    ),

    "Tissue: Glial (Astrocytes)" = list(
      genes = c("GFAP", "S100B", "AQP4", "ALDH1L1", "SLC1A3", "SLC1A2",
                "GJA1", "GLUL", "SOX9", "ALDOC")
    ),

    "Tissue: Oligodendrocytes" = list(
      genes = c("MBP", "MOG", "PLP1", "MAG", "OLIG1", "OLIG2", "SOX10",
                "CLDN11", "CNP", "MOBP")
    ),

    "Tissue: Hepatocytes" = list(
      genes = c("ALB", "APOA1", "APOB", "TTR", "FGA", "FGB", "FGG", "CYP3A4",
                "CYP2E1", "HNF4A", "SERPINA1", "TF")
    ),

    "Tissue: Kidney epithelial" = list(
      genes = c("AQP1", "AQP2", "SLC12A1", "SLC12A3", "UMOD", "NPHS1",
                "NPHS2", "PODXL", "WT1", "LRP2", "CUBN")
    ),

    "Tissue: Lung epithelial alveolar (AT1)" = list(
      genes = c("AGER", "PDPN", "HOPX", "AQP5", "CAV1", "CLIC5", "RTKN2")
    ),

    "Tissue: Lung epithelial alveolar (AT2)" = list(
      genes = c("SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "ABCA3", "LAMP3",
                "PGC", "SLC34A2", "LPCAT1")
    ),

    "Tissue: Lung epithelial ciliated" = list(
      genes = c("FOXJ1", "TUBA1A", "TUBA4A", "RSPH1", "DNAI1", "DNAH5",
                "HYDIN", "CAPS", "SNTN", "PIFO")
    ),

    "Tissue: Lung epithelial club cells" = list(
      genes = c("SCGB1A1", "SCGB3A2", "CYP2F1", "MUC5B", "CBR2", "LYPD2")
    ),

    "Tissue: Intestinal epithelial enterocytes" = list(
      genes = c("FABP1", "FABP2", "APOA1", "APOA4", "SLC2A2", "VIL1",
                "ANPEP", "DPP4", "SI", "LCT")
    ),

    "Tissue: Intestinal epithelial goblet cells" = list(
      genes = c("MUC2", "TFF3", "FCGBP", "SPINK4", "AGR2", "CLCA1", "REP15")
    ),

    "Tissue: Intestinal epithelial Paneth cells" = list(
      genes = c("LYZ", "DEFA5", "DEFA6", "PLA2G2A", "REG3A", "REG3G", "PRSS2")
    ),

    "Tissue: Intestinal stem cells" = list(
      genes = c("LGR5", "OLFM4", "ASCL2", "SMOC2", "AXIN2", "EPHB2", "SOX9")
    ),

    # ========== DEVELOPMENTAL SIGNATURES ==========

    "Development: Stem cells pluripotent" = list(
      genes = c("POU5F1", "SOX2", "NANOG", "KLF4", "MYC", "LIN28A", "DPPA4",
                "TDGF1", "ZFP42", "DNMT3B", "UTF1")
    ),

    "Development: Progenitor cells" = list(
      genes = c("CD34", "KIT", "CD133", "THY1", "ENG", "NT5E", "CD44",
                "ITGA6", "CD90")
    )
  )

  # ========================================================================
  # Enhanced Signature Library - MOUSE
  # ========================================================================

  signature_library_mouse <- list(
    "Custom" = list(genes = c()),

    # ========== IMMUNE CELL SIGNATURES - EXPANDED (MOUSE) ==========

    # B cell lineage
    "Immune: B cells" = list(
      genes = c("Cd19", "Ms4a1", "Cd79a", "Cd79b", "Pax5", "Cd22", "Bank1",
                "Blk", "Fcrl5", "Vpreb3", "Ighm", "Ighd")
    ),

    "Immune: B cells naive" = list(
      genes = c("Ighd", "Ighm", "Ccr7", "Sell", "Cd19", "Ms4a1", "Tcl1",
                "Fcer2a", "Il4r")
    ),

    "Immune: B cells memory" = list(
      genes = c("Cd27", "Cd19", "Ms4a1", "Tnfrsf13b", "Cd80", "Nt5e", "Aim2")
    ),

    "Immune: B cells activated" = list(
      genes = c("Cd69", "Cd83", "Jun", "Irf4", "Myc", "Cd86", "Tnfrsf17")
    ),

    "Immune: Plasma cells" = list(
      genes = c("Mzb1", "Sdc1", "Cd38", "Xbp1", "Prdm1", "Jchain", "Ighg1",
                "Ighg2b", "Ighg2c", "Ighg3", "Igha", "Igkc", "Iglc2", "Iglc3",
                "Ssr4", "Derl3", "Fkbp11", "Prdx4", "Slamf7", "Cd27")
    ),

    "Immune: Plasmablasts" = list(
      genes = c("Ighg1", "Ighg2b", "Jchain", "Iglc2", "Mzb1", "Igkc", "Ighm",
                "Xbp1", "Cd38", "Cd27", "Sdc1")
    ),

    # T cell lineage
    "Immune: T cells CD8" = list(
      genes = c("Cd8a", "Cd8b1", "Cd3d", "Cd3e", "Cd3g", "Gzmk", "Ccl5",
                "Gzma", "Prf1", "Nkg7", "Trac", "Trbc1", "Trbc2")
    ),

    "Immune: T cells CD8 naive" = list(
      genes = c("Cd8a", "Cd8b1", "Ccr7", "Sell", "Lef1", "Tcf7", "Cd27",
                "Cd28", "Il7r")
    ),

    "Immune: T cells CD8 effector" = list(
      genes = c("Cd8a", "Gzmb", "Gzmh", "Prf1", "Gnly", "Nkg7", "Fgfbp2",
                "Fcgr3", "Cx3cr1")
    ),

    "Immune: T cells CD8 exhausted" = list(
      genes = c("Cd8a", "Pdcd1", "Ctla4", "Havcr2", "Lag3", "Tigit", "Tox",
                "Entpd1", "Layn", "Cxcl13")
    ),

    "Immune: T cells CD4" = list(
      genes = c("Cd4", "Il7r", "Cd3d", "Cd3e", "Cd3g", "Cd40lg", "Icos",
                "Maf", "Trac", "Trbc1", "Trbc2", "Cd28")
    ),

    "Immune: T cells CD4 naive" = list(
      genes = c("Cd4", "Ccr7", "Sell", "Lef1", "Tcf7", "Cd27", "Cd28",
                "Il7r", "Mal")
    ),

    "Immune: T cells CD4 memory" = list(
      genes = c("Cd4", "Il7r", "Ccr7", "Sell", "Cd27", "Cd28", "Tnfrsf4",
                "Ltb", "Aqp3")
    ),

    "Immune: T cells regulatory (Tregs)" = list(
      genes = c("Foxp3", "Il2ra", "Ctla4", "Ikzf2", "Tnfrsf18", "Tnfrsf4",
                "Il10", "Tigit", "Layn", "Ccr8", "Il2rb", "Cd4")
    ),

    "Immune: T cells follicular helper (Tfh)" = list(
      genes = c("Cxcr5", "Pdcd1", "Bcl6", "Cd4", "Il21", "Icos", "Cd40lg",
                "Maf", "Tox2", "Btla")
    ),

    "Immune: T cells Th1" = list(
      genes = c("Cd4", "Ifng", "Tbx21", "Il12rb2", "Cxcr3", "Ccr5", "Tnf")
    ),

    "Immune: T cells Th2" = list(
      genes = c("Cd4", "Gata3", "Il4", "Il5", "Il13", "Ccr4", "Il1rl1")
    ),

    "Immune: T cells Th17" = list(
      genes = c("Cd4", "Rorc", "Il17a", "Il17f", "Il23r", "Ccr6", "Il22")
    ),

    # NK cells
    "Immune: NK cells" = list(
      genes = c("Ncam1", "Fcgr3", "Nkg7", "Gzma", "Klrd1", "Klrf1", "Klrb1c",
                "Klrc1", "Klrc2", "Klrc3", "Ncr1", "Gzmb", "Prf1", "Xcl1",
                "Fcer1g", "Tyrobp")
    ),

    "Immune: NK cells activated" = list(
      genes = c("Ncam1", "Sell", "Il7r", "Cd44", "Xcl1", "Gzmk")
    ),

    "Immune: NK cells cytotoxic" = list(
      genes = c("Fcgr3", "Klrf1", "Klrd1", "Prf1", "Gzmb", "S1pr5", "Cx3cr1")
    ),

    # Myeloid lineage
    "Immune: Monocytes" = list(
      genes = c("Cd14", "Fcgr3", "S100a8", "S100a9", "Lyz2", "Vcan", "Fcn1",
                "S100a4", "Csf1r", "Cd68", "Itgam", "Itgax")
    ),

    "Immune: Monocytes classical" = list(
      genes = c("Ly6c2", "Ccr2", "S100a8", "S100a9", "S100a4", "Lyz2",
                "Fcn1", "Vcan", "Itgam")
    ),

    "Immune: Monocytes non-classical" = list(
      genes = c("Fcgr3", "Cd14", "Cdkn1c", "Ms4a7", "Asah1", "Ifitm3",
                "Lst1", "Nr4a1")
    ),

    "Immune: Macrophages" = list(
      genes = c("Cd68", "Cd163", "Csf1r", "Msr1", "C1qa", "C1qb", "C1qc",
                "Marco", "Mrc1", "Cd14", "Fcgr1", "Fcgr2b")
    ),

    "Immune: Macrophages M1" = list(
      genes = c("Cd68", "Cd80", "Cd86", "Il1b", "Il6", "Tnf", "Cxcl9",
                "Cxcl10", "Nos2", "Ido1", "Socs3")
    ),

    "Immune: Macrophages M2" = list(
      genes = c("Cd68", "Cd163", "Mrc1", "Msr1", "Lyve1", "Maf", "Folr2",
                "F13a1", "Sepp1", "Stab1", "Cd209a", "Clec10a")
    ),

    "Immune: Macrophages alveolar" = list(
      genes = c("Marco", "Fabp4", "Pparg", "Mcemp1", "Cd68", "Slc7a8",
                "Chit1", "Msr1", "Car4")
    ),

    "Immune: Macrophages tumor-associated (TAM)" = list(
      genes = c("Cd68", "Cd163", "Mrc1", "Msr1", "Trem2", "Spp1", "Apoe",
                "C1qa", "C1qb", "Folr2", "Lyve1")
    ),

    "Immune: Dendritic cells" = list(
      genes = c("H2-Aa", "H2-Ab1", "H2-Eb1", "Cd1c", "Fcer1a", "Clec10a",
                "Cd209a", "Itgax", "Itgam")
    ),

    "Immune: Dendritic cells conventional (cDC1)" = list(
      genes = c("Clec9a", "Xcr1", "Cadm1", "Batf3", "Irf8", "Ido1", "Itgae")
    ),

    "Immune: Dendritic cells conventional (cDC2)" = list(
      genes = c("Cd1c", "Fcer1a", "Clec10a", "Cd209a", "Sirpa", "Irf4")
    ),

    "Immune: Dendritic cells plasmacytoid (pDC)" = list(
      genes = c("Il3ra", "Clec4c", "Nrp1", "Tcf4", "Irf7", "Irf8", "Ly6d",
                "Gzmb", "Itm2c", "Siglech")
    ),

    "Immune: Dendritic cells activated" = list(
      genes = c("Cd83", "Lamp3", "Ccr7", "Fscn1", "Ccl19", "Ccl22", "Il1b",
                "Marcksl1", "Ido1")
    ),

    "Immune: Neutrophils" = list(
      genes = c("Ly6g", "Csf3r", "S100a8", "S100a9", "S100a4", "Cxcr2",
                "Fpr1", "Fcgr2b", "Mpo", "Elane", "Ctsg", "Prtn3", "Lcn2",
                "Mmp8", "Mmp9", "Camp", "Ltf")
    ),

    "Immune: Eosinophils" = list(
      genes = c("Siglecf", "Il5ra", "Ccr3", "Ear2", "Prg2", "Prg3",
                "Rnase2a", "Alox15", "Il1rl1", "Epx")
    ),

    "Immune: Basophils" = list(
      genes = c("Mcpt8", "Hdc", "Ms4a2", "Gata2", "Cpa3", "Tpsab1", "Il3ra",
                "Slc18a2", "Il4", "Il13")
    ),

    "Immune: Mast cells" = list(
      genes = c("Tpsab1", "Tpsb2", "Cpa3", "Kit", "Gata2", "Hdc", "Ms4a2",
                "Hpgds", "Ctsg", "Rgs13", "Il1rl1", "Fcer1a", "Slc18a2", "Mcpt4")
    ),

    # ========== FUNCTIONAL IMMUNE SIGNATURES - EXPANDED (MOUSE) ==========

    "Immune: Cytotoxicity" = list(
      genes = c("Prf1", "Gzma", "Gzmb", "Gzmh", "Gzmk", "Gzmc", "Gnly",
                "Nkg7", "Klrk1", "Klrb1c", "Fasl", "Ifng", "Tnf")
    ),

    "Immune: Exhaustion" = list(
      genes = c("Pdcd1", "Ctla4", "Havcr2", "Lag3", "Tigit", "Tox", "Tox2",
                "Entpd1", "Itgae", "Layn", "Cxcl13", "Cd244", "Cd160",
                "Btla", "Vsir", "Tnfrsf9")
    ),

    "Immune: Inflammation" = list(
      genes = c("Il1b", "Il6", "Tnf", "Cxcl1", "Cxcl2", "Il1a", "Ptgs2",
                "Ccl2", "Ccl3", "Ccl4", "Cxcl3", "Csf2", "Csf3")
    ),

    "Immune: Interferon response" = list(
      genes = c("Ifitm1", "Ifitm2", "Ifitm3", "Isg15", "Isg20", "Ifi27l2a",
                "Ifi44", "Ifi44l", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2",
                "Oas1a", "Oas2", "Oas3", "Irf7", "Stat1", "Stat2")
    ),

    "Immune: Complement activation" = list(
      genes = c("C1qa", "C1qb", "C1qc", "C1r", "C1s", "C2", "C3", "C4a",
                "C4b", "C5", "Cfb", "Cfd", "Cfh", "Cfi")
    ),

    "Immune: Phagocytosis" = list(
      genes = c("Fcgr1", "Fcgr2b", "Fcgr3", "Fcgr4", "Msr1", "Cd36", "Marco",
                "Scarb1", "Mertk", "Axl", "Tyro3")
    ),

    "Immune: Antigen presentation" = list(
      genes = c("H2-K1", "H2-D1", "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74",
                "B2m", "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9")
    ),

    "Immune: Chemotaxis" = list(
      genes = c("Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5", "Ccr6", "Ccr7",
                "Ccr9", "Cxcr1", "Cxcr2", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6",
                "Cx3cr1", "Ackr1", "Ackr2")
    ),

    # ========== CANCER SIGNATURES - EXPANDED (MOUSE) ==========

    "Cancer: Epithelial" = list(
      genes = c("Epcam", "Krt8", "Krt18", "Krt19", "Cdh1", "Cldn3", "Cldn4",
                "Cldn7", "Tacstd2", "Muc1", "Elf3", "Krt7", "Krt20", "Grhl2")
    ),

    "Cancer: Proliferation" = list(
      genes = c("Mki67", "Top2a", "Pcna", "Cdk1", "Ccnb1", "Ccnb2", "Ccne1",
                "Ccne2", "Cdc20", "Aurka", "Aurkb", "Plk1", "Bub1", "Bub1b",
                "Cenpe", "Cenpf", "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6",
                "Mcm7", "Tyms", "Tk1", "Rrm2", "Pttg1", "Ube2c", "Birc5")
    ),

    "Cancer: Hypoxia" = list(
      genes = c("Hif1a", "Vegfa", "Car9", "Slc2a1", "Ldha", "Pgk1", "Bnip3",
                "Bnip3l", "Pdk1", "Eno1", "Pfkl", "Pfkp", "Slc16a3", "Ndrg1",
                "Egln3", "Lox", "P4ha1", "P4ha2", "Adm", "Angptl4", "Vhl")
    ),

    "Cancer: EMT (Epithelial-Mesenchymal Transition)" = list(
      genes = c("Vim", "Cdh2", "Snai1", "Snai2", "Twist1", "Twist2", "Zeb1",
                "Zeb2", "Fn1", "Mmp2", "Mmp3", "Mmp9", "Mmp14", "Tcf4",
                "Foxc2", "Sparc", "Cdh11", "Thy1", "Itga5", "Itgb1")
    ),

    "Cancer: Stemness" = list(
      genes = c("Cd44", "Aldh1a1", "Prom1", "Sox2", "Nanog", "Pou5f1",
                "Klf4", "Myc", "Lgr5", "Bmi1", "Notch1", "Notch2", "Jag1",
                "Dll1", "Cd24a", "Epcam", "Itga6", "Thy1")
    ),

    "Cancer: Angiogenesis" = list(
      genes = c("Vegfa", "Vegfb", "Vegfc", "Kdr", "Flt1", "Flt4", "Pecam1",
                "Angpt1", "Angpt2", "Tek", "Nrp1", "Nrp2", "Pgf", "Hif1a",
                "Pdgfa", "Pdgfb", "Pdgfra", "Pdgfrb", "Fgf2", "Fgfr1")
    ),

    "Cancer: Metastasis" = list(
      genes = c("Mmp2", "Mmp9", "Mmp11", "Mmp14", "Cxcr4", "Cxcl12", "Twist1",
                "Snai1", "Zeb1", "Vim", "Fn1", "Cdh2", "Itga5", "Itgb3",
                "Itgav", "Cd44", "Met", "Ptgs2", "Ctnnb1")
    ),

    "Cancer: Apoptosis resistance" = list(
      genes = c("Bcl2", "Bcl2l1", "Mcl1", "Birc5", "Xiap", "Birc2", "Birc3",
                "Mdm2", "Trp53", "Bax", "Bak1", "Bid", "Casp3", "Casp8",
                "Casp9", "Cycs")
    ),

    "Cancer: DNA damage response" = list(
      genes = c("Brca1", "Brca2", "Atm", "Atr", "Chek1", "Chek2", "Rad51",
                "Xrcc1", "Parp1", "H2afx", "Tp53bp1", "Mdc1", "Fancd2",
                "Fanca", "Fancg")
    ),

    # ========== STROMAL SIGNATURES - EXPANDED (MOUSE) ==========

    "Stroma: Fibroblasts" = list(
      genes = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra", "Thy1",
                "Vim", "Acta2", "Tagln", "Myl9", "Fn1", "Sparc")
    ),

    "Stroma: Fibroblasts (CAF)" = list(
      genes = c("Fap", "Pdgfra", "Pdpn", "Col1a1", "Col1a2", "Col3a1",
                "Col5a1", "Col6a1", "Fn1", "Dcn", "Lum", "Acta2", "Tagln",
                "Tgfb1", "Tgfb2", "Cxcl12", "Cxcl14", "Il6", "Ccl2",
                "Postn", "Thbs2", "Mmp11", "Lrrc15", "Comp")
    ),

    "Stroma: Myofibroblasts" = list(
      genes = c("Acta2", "Tagln", "Myl9", "Myh11", "Cnn1", "Actg2", "Des",
                "Pdgfrb", "Rgs5")
    ),

    "Stroma: Endothelial" = list(
      genes = c("Pecam1", "Vwf", "Cdh5", "Kdr", "Flt1", "Eng", "Mcam",
                "Plvap", "Tek", "Esam", "Cldn5", "Cd34", "Erg", "Fli1",
                "Icam2", "Emcn", "Nos3", "Ramp2", "Ramp3")
    ),

    "Stroma: Endothelial arterial" = list(
      genes = c("Gja5", "Bmx", "Efnb2", "Sox17", "Dll4", "Hey1", "Sema3g")
    ),

    "Stroma: Endothelial venous" = list(
      genes = c("Nr2f2", "Ackr1", "Selp", "Vwf", "Clu", "Aplnr")
    ),

    "Stroma: Endothelial lymphatic" = list(
      genes = c("Prox1", "Lyve1", "Pdpn", "Flt4", "Ccl21a", "Tff3", "Mmrn1")
    ),

    "Stroma: Endothelial tip cells" = list(
      genes = c("Apln", "Esm1", "Angpt2", "Flt1", "Apln", "Pgf", "Kdr",
                "Pdgfb", "Dll4")
    ),

    "Stroma: Pericytes" = list(
      genes = c("Pdgfrb", "Rgs5", "Acta2", "Mcam", "Notch3", "Cspg4",
                "Abcc9", "Kcnj8", "Ndufa4l2", "Cox4i2", "Steap4", "Angpt1")
    ),

    "Stroma: Smooth muscle cells" = list(
      genes = c("Acta2", "Myh11", "Tagln", "Cnn1", "Actg2", "Myl9", "Des",
                "Mylk", "Smtn", "Tpm2", "Tpm1", "Lmod1")
    ),

    "Stroma: Adipocytes" = list(
      genes = c("Adipoq", "Lep", "Plin1", "Fabp4", "Lpl", "Pparg", "Cebpa",
                "Cidec", "Plin4", "Lipe", "Dgat2", "Slc2a4")
    ),

    # ========== TISSUE-SPECIFIC SIGNATURES (MOUSE) ==========

    "Tissue: Neural" = list(
      genes = c("Tubb3", "Map2", "Syp", "Snap25", "Syt1", "Nefl", "Eno2",
                "Dcx", "Ncam1", "Gap43", "Stmn2", "Rbfox3")
    ),

    "Tissue: Glial (Astrocytes)" = list(
      genes = c("Gfap", "S100b", "Aqp4", "Aldh1l1", "Slc1a3", "Slc1a2",
                "Gja1", "Glul", "Sox9", "Aldoc")
    ),

    "Tissue: Oligodendrocytes" = list(
      genes = c("Mbp", "Mog", "Plp1", "Mag", "Olig1", "Olig2", "Sox10",
                "Cldn11", "Cnp", "Mobp")
    ),

    "Tissue: Hepatocytes" = list(
      genes = c("Alb", "Apoa1", "Apob", "Ttr", "Fga", "Fgb", "Fgg", "Cyp3a11",
                "Cyp2e1", "Hnf4a", "Serpina1a", "Tf")
    ),

    "Tissue: Kidney epithelial" = list(
      genes = c("Aqp1", "Aqp2", "Slc12a1", "Slc12a3", "Umod", "Nphs1",
                "Nphs2", "Podxl", "Wt1", "Lrp2", "Cubn")
    ),

    "Tissue: Lung epithelial alveolar (AT1)" = list(
      genes = c("Ager", "Pdpn", "Hopx", "Aqp5", "Cav1", "Clic5", "Rtkn2")
    ),

    "Tissue: Lung epithelial alveolar (AT2)" = list(
      genes = c("Sftpc", "Sftpb", "Sftpa1", "Sftpd", "Abca3", "Lamp3",
                "Pgc", "Slc34a2", "Lpcat1")
    ),

    "Tissue: Lung epithelial ciliated" = list(
      genes = c("Foxj1", "Tuba1a", "Tuba1b", "Rsph1", "Dnai1", "Dnah5",
                "Hydin", "Caps", "Sntn", "Pifo")
    ),

    "Tissue: Lung epithelial club cells" = list(
      genes = c("Scgb1a1", "Scgb3a2", "Cyp2f2", "Muc5b", "Cbr2", "Lypd2")
    ),

    "Tissue: Intestinal epithelial enterocytes" = list(
      genes = c("Fabp1", "Fabp2", "Apoa1", "Apoa4", "Slc2a2", "Vil1",
                "Anpep", "Dpp4", "Si", "Lct")
    ),

    "Tissue: Intestinal epithelial goblet cells" = list(
      genes = c("Muc2", "Tff3", "Fcgbp", "Spink4", "Agr2", "Clca1", "Zg16")
    ),

    "Tissue: Intestinal epithelial Paneth cells" = list(
      genes = c("Lyz1", "Defa24", "Defa-rs1", "Pla2g2a", "Ang4", "Prss2")
    ),

    "Tissue: Intestinal stem cells" = list(
      genes = c("Lgr5", "Olfm4", "Ascl2", "Smoc2", "Axin2", "Ephb2", "Sox9")
    ),

    # ========== DEVELOPMENTAL SIGNATURES (MOUSE) ==========

    "Development: Stem cells pluripotent" = list(
      genes = c("Pou5f1", "Sox2", "Nanog", "Klf4", "Myc", "Lin28a", "Dppa4",
                "Tdgf1", "Zfp42", "Dnmt3b", "Utf1")
    ),

    "Development: Progenitor cells" = list(
      genes = c("Cd34", "Kit", "Prom1", "Thy1", "Eng", "Nt5e", "Cd44",
                "Itga6", "Cd90")
    )
  )

  # Return everything
  list(
    seurat_obj = seurat_obj,
    coords = coords,
    image_name = image_name,
    spots_sf = spots_sf,
    he_image_base64 = he_image_base64,
    he_image_bounds = he_image_bounds,
    signature_library_human = signature_library_human,
    signature_library_mouse = signature_library_mouse,
    sample_name = sample_name,
    show_image = show_image
  )
}
