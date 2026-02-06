#!/usr/bin/env Rscript

# Define function to move a column in a data table
move_dt_column <- function(table, column_being_moved, column_before) {
    setcolorder(table, c(
        names(table)[1:which(names(table) == column_before)],
        column_being_moved,
        setdiff(names(table), c(names(table)[1:which(names(table) == column_before)], column_being_moved))
    ))
    return(table)
}

# Define function to clean variant data table
clean_variant_table <- function(table, ancestry, cadd_deleterious_threshold, conservation_conserved_threshold, kinase_activation_threshold, modality) {
    table <- as.data.table(table)

    # Compute frequencies
    if (all(c("het_PD", "hom_PD", "total_PD") %in% names(table))) {
        table[, `PD Frequency` := (het_PD + 2*hom_PD) / (2*total_PD)]
    }
    if (all(c("het_HC", "hom_HC", "total_HC") %in% names(table))) {
        table[, `Control Frequency` := (het_HC + 2*hom_HC) / (2*total_HC)]
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
        SAS = "gnomad41_genome_AF_sas",
        AAC = NA,
        CAS = NA,
        CAH = NA,
        Combined = "gnomad41_genome_AF"
    )

    # Use general gnomad frequency for missing ancestries
    gnomad_col <- ancestry_col_map[[ancestry]]
    if (is.na(gnomad_col)) {
        table[, "gnomAD allele frequency" := NA_real_]
        gnomad_col <- "gnomAD allele frequency"
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
        "Clinvar_Pathogenic",
        "exon",
        "cDNA",
        "prot_change",
        "domain",
        "Consurf_score",
        "Mean_pRAB10/RAB10",
        "FH_neg_PD",
        "FH_pos_PD",
        "FH_unk_PD",
        "FH_neg_HC",
        "FH_pos_HC",
        "FH_unk_HC",
        "Age_10",
        "Age_20",
        "Age_30",
        "Age_40",
        "Age_50",
        "Age_60",
        "Age_70",
        "Age_80",
        "Age_90",
        "Age_100",
        "Age_min",
        "Age_max",
        "Age_median"
    )
    if (modality == "Clinical exome") {
        keep_cols <- keep_cols[!keep_cols %in% c("Control Frequency")]
    }

    table <- table[, intersect(keep_cols, names(table)), with = FALSE]

    # Replace "." with NA
    na_cols <- c("ExonicFunc.refGene", "CADD_phred", "eQTLGen_snp_id", gnomad_col, "CLNSIG", "CLNDN")
    for (col in intersect(na_cols, names(table))) {
        table[[col]][table[[col]] == "."] <- NA
    }

    # Convert numeric columns
    for (col in c("CADD_phred", gnomad_col)) {
        if (col %in% names(table)) {
            table[[col]] <- suppressWarnings(as.numeric(table[[col]]))
        }
    }

    # Clean up exon column
    if ("exon" %in% names(table)) {
        table$exon <- gsub("exon", "", table$exon)
        table$exon <- suppressWarnings(as.numeric(table$exon))
    }

    # Clean up clinical significance
    table[, ("Clinvar_Pathogenic") := lapply(.SD, function(x) {
        x <- gsub("_", " ", x, fixed = TRUE)
        x <- gsub("|", ", ", x, fixed = TRUE)
        x
    }), .SDcols = "Clinvar_Pathogenic"]

    # Rename columns
    new_col_names <- c(
        "Variant (GrCh38)",
        paste0("PD frequency (", modality, ")"),
        paste0("Control frequency (", modality, ")"),
        "gnomAD allele frequency",
        "Region",
        "Functional consequence",
        "CADD",
        "rsID",
        "Clinvar Pathogenic",
        "Exon #",
        "cDNA change",
        "AA change",
        "Protein domain",
        "Conservation score",
        "Kinase activity (mean pRAB10/RAB10)",
        "No family history (PD)",
        "Family history (PD)",
        "Unknown family history (PD)",
        "No family history (Control)",
        "Family history (Control)",
        "Unknown family history (Control)",
        "AAO (1-10)",
        "AAO (11-20)",
        "AAO (21-30)",
        "AAO (31-40)",
        "AAO (41-50)",
        "AAO (51-60)",
        "AAO (61-70)",
        "AAO (71-80)",
        "AAO (81-90)",
        "AAO (91-100)",
        "Minimum AAO",
        "Maximum AAO",
        "Median AAO"
    )
    if (modality == "Clinical exome") {
        new_col_names <- new_col_names[!new_col_names %in% c(paste0("Control frequency (", modality, ")"))]
    }

    colnames(table) <- new_col_names

    table[, `Deleterious?` := fifelse(CADD > cadd_deleterious_threshold, "Yes", "")]
    table[, `Conserved?` := fifelse(`Conservation score` > conservation_conserved_threshold, "Yes", "")]
    table[, `Kinase active?` := fifelse(`Kinase activity (mean pRAB10/RAB10)` > kinase_activation_threshold, "Yes", "")]

    table <- move_dt_column(table, "Deleterious?", "CADD")
    table <- move_dt_column(table, "Conserved?", "Conservation score")
    table <- move_dt_column(table, "Kinase active?", "Kinase activity (mean pRAB10/RAB10)")

    return(table)
}
