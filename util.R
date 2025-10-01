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

# Define function to process tables
clean_variant_table <- function(tbl, cadd_deleterious_threshold, conservation_conserved_threshold, kinase_activation_threshold) {
    # Find frequency of cases and controls
    tbl <- as.data.table(tbl)
    if (all(c("het_PD", "hom_PD", "total_PD") %in% names(tbl))) {
        tbl[, `PD Frequency` := (het_PD + 2*hom_PD) / (2*total_PD)]
    }
    if (all(c("het_HC", "hom_HC", "total_HC") %in% names(tbl))) {
        tbl[, `Control Frequency` := (het_HC + 2*hom_HC) / (2*total_HC)]
    }

    # Keep only selected columns
    keep_cols <- c(
        "variant",
        "PD Frequency",
        "Control Frequency",
        "gnomad41_genome_AF",
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
    )
    tbl <- tbl[, intersect(keep_cols, names(tbl)), with = FALSE]

    # Replace "." with NA so they appear blank in the table
    na_cols <- c("ExonicFunc.refGene", "CADD_phred", "eQTLGen_snp_id", "gnomad41_genome_AF", "CLNSIG", "CLNDN")
    for (col in intersect(na_cols, names(tbl))) {
        tbl[[col]][tbl[[col]] == "."] <- NA
    }

    # Convert numeric columns
    for (col in c("CADD_phred", "gnomad41_genome_AF")) {
        if (col %in% names(tbl)) {
            tbl[[col]] <- suppressWarnings(as.numeric(tbl[[col]]))
        }
    }

    # Clean up exon column
    if ("exon" %in% names(tbl)) {
        tbl$exon <- gsub("exon", "", tbl$exon)
        tbl$exon <- suppressWarnings(as.numeric(tbl$exon))
    }

    # Clean up clinical significance columns
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
    )

    tbl[, `Deleterious?` := fifelse(CADD > cadd_deleterious_threshold, "Yes", "")]
    tbl[, `Conserved?` := fifelse(`Conservation score` > conservation_conserved_threshold, "Yes", "")]
    tbl[, `Kinase active?` := fifelse(`Kinase activity (mean pRAB10/RAB10)` > kinase_activation_threshold, "Yes", "")]

    tbl <- move_dt_column(tbl, "Deleterious?", "CADD")
    tbl <- move_dt_column(tbl, "Conserved?", "Conservation score")
    tbl <- move_dt_column(tbl, "Kinase active?", "Kinase activity (mean pRAB10/RAB10)")

    return(tbl)
}
