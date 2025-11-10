#!/usr/bin/env Rscript

# Define function to move a column in a data table
move_dt_column <- function(tbl, column_being_moved, column_before) {
    setcolorder(tbl, c(
        names(tbl)[1:which(names(tbl) == column_before)],
        column_being_moved,
        setdiff(names(tbl), c(names(tbl)[1:which(names(tbl) == column_before)], column_being_moved))
    ))
    return(tbl)
}

# # Define function to process tables
# clean_variant_table <- function(tbl, cadd_deleterious_threshold, conservation_conserved_threshold, kinase_activation_threshold) {
#     # Find frequency of cases and controls
#     tbl <- as.data.table(tbl)
#     if (all(c("het_PD", "hom_PD", "total_PD") %in% names(tbl))) {
#         tbl[, `PD Frequency` := (het_PD + 2*hom_PD) / (2*total_PD)]
#     }
#     if (all(c("het_HC", "hom_HC", "total_HC") %in% names(tbl))) {
#         tbl[, `Control Frequency` := (het_HC + 2*hom_HC) / (2*total_HC)]
#     }

#     # Keep only selected columns
#     keep_cols <- c(
#         "variant",
#         "PD Frequency",
#         "Control Frequency",
#         "gnomad41_genome_AF",
#         "gnomad41_genome_AF_afr",
#         "gnomad41_genome_AF_amr",
#         "gnomad41_genome_AF_asj",
#         "gnomad41_genome_AF_eas",
#         "gnomad41_genome_AF_fin",
#         "gnomad41_genome_AF_mid",
#         "gnomad41_genome_AF_nfe",
#         "gnomad41_genome_AF_remaining",
#         "gnomad41_genome_AF_sas",
#         "Func.refGene",
#         "ExonicFunc.refGene",
#         "CADD_phred",
#         "eQTLGen_snp_id",
#         "CLNSIG",
#         "CLNDN",
#         "exon",
#         "cDNA",
#         "prot_change",
#         "domain",
#         "Consurf_score",
#         "Mean_pRAB10/RAB10",
#         "FH_neg_PD",
#         "FH_pos_PD",
#         "FH_unk_PD",
#         "FH_neg_HC",
#         "FH_pos_HC",
#         "FH_unk_HC",
#         "Age_10",
#         "Age_20",
#         "Age_30",
#         "Age_40",
#         "Age_50",
#         "Age_60",
#         "Age_70",
#         "Age_80",
#         "Age_90",
#         "Age_100",
#         "Age_min",
#         "Age_max",
#         "Age_median"
#     )
#     tbl <- tbl[, intersect(keep_cols, names(tbl)), with = FALSE]

#     # Replace "." with NA so they appear blank in the table
#     na_cols <- c("ExonicFunc.refGene", "CADD_phred", "eQTLGen_snp_id", "gnomad41_genome_AF", "CLNSIG", "CLNDN")
#     for (col in intersect(na_cols, names(tbl))) {
#         tbl[[col]][tbl[[col]] == "."] <- NA
#     }

#     # Convert numeric columns
#     for (col in c("CADD_phred", "gnomad41_genome_AF")) {
#         if (col %in% names(tbl)) {
#             tbl[[col]] <- suppressWarnings(as.numeric(tbl[[col]]))
#         }
#     }

#     # Clean up exon column
#     if ("exon" %in% names(tbl)) {
#         tbl$exon <- gsub("exon", "", tbl$exon)
#         tbl$exon <- suppressWarnings(as.numeric(tbl$exon))
#     }

#     # Clean up clinical significance columns
#     tbl[, (c("CLNSIG","CLNDN")) := lapply(.SD, function(x) {
#         x <- gsub("_", " ", x, fixed = TRUE)
#         x <- gsub("|", ", ", x, fixed = TRUE)
#         x
#     }), .SDcols = c("CLNSIG","CLNDN")]

#     # Rename columns
#     colnames(tbl) <- c(
#         "Variant (GrCh38)",
#         "PD frequency",
#         "Control frequency",
#         "Gnomad allele frequency (Combined)",
#         "Gnomad allele frequency (AFR)",
#         "Gnomad allele frequency (AMR)",
#         "Gnomad allele frequency (AJ)",
#         "Gnomad allele frequency (EAS)",
#         "Gnomad allele frequency (FIN)",
#         "Gnomad allele frequency (MDE)",
#         "Gnomad allele frequency (EUR)",
#         "Gnomad allele frequency (CAH)",
#         "Gnomad allele frequency (SAS)",
#         "Region",
#         "Functional consequence",
#         "CADD",
#         "rsID",
#         "Clinical significance",
#         "Clinical disease name",
#         "Exon #",
#         "cDNA change",
#         "AA change",
#         "Protein domain",
#         "Conservation score",
#         "Kinase activity (mean pRAB10/RAB10)",
#         "FH_neg_PD",
#         "FH_pos_PD",
#         "FH_unk_PD",
#         "FH_neg_HC",
#         "FH_pos_HC",
#         "FH_unk_HC",
#         "Age_10",
#         "Age_20",
#         "Age_30",
#         "Age_40",
#         "Age_50",
#         "Age_60",
#         "Age_70",
#         "Age_80",
#         "Age_90",
#         "Age_100",
#         "Age_min",
#         "Age_max",
#         "Age_median"
#     )

#     tbl[, `Deleterious?` := fifelse(CADD > cadd_deleterious_threshold, "Yes", "")]
#     tbl[, `Conserved?` := fifelse(`Conservation score` > conservation_conserved_threshold, "Yes", "")]
#     tbl[, `Kinase active?` := fifelse(`Kinase activity (mean pRAB10/RAB10)` > kinase_activation_threshold, "Yes", "")]

#     tbl <- move_dt_column(tbl, "Deleterious?", "CADD")
#     tbl <- move_dt_column(tbl, "Conserved?", "Conservation score")
#     tbl <- move_dt_column(tbl, "Kinase active?", "Kinase activity (mean pRAB10/RAB10)")

#     return(tbl)
# }


clean_variant_table <- function(tbl, ancestry, cadd_deleterious_threshold, conservation_conserved_threshold, kinase_activation_threshold) {
    tbl <- as.data.table(tbl)

    # Compute frequencies
    if (all(c("het_PD", "hom_PD", "total_PD") %in% names(tbl))) {
        tbl[, `PD Frequency` := (het_PD + 2*hom_PD) / (2*total_PD)]
    }
    if (all(c("het_HC", "hom_HC", "total_HC") %in% names(tbl))) {
        tbl[, `Control Frequency` := (het_HC + 2*hom_HC) / (2*total_HC)]
    }

    # Identify correct gnomAD column to use
    ancestry_col_map <- list(
        AFR = "gnomad41_genome_AF_afr",
        AMR = "gnomad41_genome_AF_amr",
        AJ = "gnomad41_genome_AF_asj",
        EAS = "gnomad41_genome_AF_eas",
        FIN = "gnomad41_genome_AF_fin",
        MDE = "gnomad41_genome_AF_mid",
        EUR = "gnomad41_genome_AF_nfe",
        CAH = "gnomad41_genome_AF_remaining",
        SAS = "gnomad41_genome_AF_sas",
        Combined = "gnomad41_genome_AF"
    )

    # Use general gnomad frequency for missing ancestries
    gnomad_col <- ancestry_col_map[[ancestry]]
    if (is.null(gnomad_col) || !(gnomad_col %in% names(tbl))) {
        gnomad_col <- "gnomad41_genome_AF"
    }

    # Keep only relevant columns (including chosen gnomad col)
    keep_cols <- c(
        "variant",
        "PD Frequency",
        "Control Frequency",
        gnomad_col,
        "Func.refGene",
        "ExonicFunc.refGene",
        "CADD_phred",
        "eQTLGen_snp_id",
        "CLNSIG",
        "CLNDN",
        "exon",
        "cDNA",
        "prot_change",
        "domain",
        "Consurf_score",
        "Mean_pRAB10/RAB10"
        # "FH_neg_PD",
        # "FH_pos_PD",
        # "FH_unk_PD",
        # "FH_neg_HC",
        # "FH_pos_HC",
        # "FH_unk_HC",
        # "Age_10",
        # "Age_20",
        # "Age_30",
        # "Age_40",
        # "Age_50",
        # "Age_60",
        # "Age_70",
        # "Age_80",
        # "Age_90",
        # "Age_100",
        # "Age_min",
        # "Age_max",
        # "Age_median"
    )
    tbl <- tbl[, intersect(keep_cols, names(tbl)), with = FALSE]

    # Replace "." with NA
    na_cols <- c("ExonicFunc.refGene", "CADD_phred", "eQTLGen_snp_id", gnomad_col, "CLNSIG", "CLNDN")
    for (col in intersect(na_cols, names(tbl))) {
        tbl[[col]][tbl[[col]] == "."] <- NA
    }

    # Convert numeric columns
    for (col in c("CADD_phred", gnomad_col)) {
        if (col %in% names(tbl)) {
            tbl[[col]] <- suppressWarnings(as.numeric(tbl[[col]]))
        }
    }

    # Clean up exon column
    if ("exon" %in% names(tbl)) {
        tbl$exon <- gsub("exon", "", tbl$exon)
        tbl$exon <- suppressWarnings(as.numeric(tbl$exon))
    }

    # Clean up clinical significance
    tbl[, (c("CLNSIG","CLNDN")) := lapply(.SD, function(x) {
        x <- gsub("_", " ", x, fixed = TRUE)
        x <- gsub("|", ", ", x, fixed = TRUE)
        x
    }), .SDcols = c("CLNSIG","CLNDN")]

    # Rename columns
    colnames(tbl) <- c(
        "Variant (GrCh38)",
        "PD frequency",
        "Control frequency",
        "Gnomad allele frequency",
        "Region",
        "Functional consequence",
        "CADD",
        "rsID",
        "Clinical significance",
        "Clinical disease name",
        "Exon #",
        "cDNA change",
        "AA change",
        "Protein domain",
        "Conservation score",
        "Kinase activity (mean pRAB10/RAB10)"
        # "FH_neg_PD",
        # "FH_pos_PD",
        # "FH_unk_PD",
        # "FH_neg_HC",
        # "FH_pos_HC",
        # "FH_unk_HC",
        # "Age_10",
        # "Age_20",
        # "Age_30",
        # "Age_40",
        # "Age_50",
        # "Age_60",
        # "Age_70",
        # "Age_80",
        # "Age_90",
        # "Age_100",
        # "Age_min",
        # "Age_max",
        # "Age_median"
    )

    tbl[, `Deleterious?` := fifelse(CADD > cadd_deleterious_threshold, "Yes", "")]
    tbl[, `Conserved?` := fifelse(`Conservation score` > conservation_conserved_threshold, "Yes", "")]
    tbl[, `Kinase active?` := fifelse(`Kinase activity (mean pRAB10/RAB10)` > kinase_activation_threshold, "Yes", "")]

    tbl <- move_dt_column(tbl, "Deleterious?", "CADD")
    tbl <- move_dt_column(tbl, "Conserved?", "Conservation score")
    tbl <- move_dt_column(tbl, "Kinase active?", "Kinase activity (mean pRAB10/RAB10)")

    return(tbl)
}
