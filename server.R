#!/usr/bin/env Rscript

# =========================================================================
# INITIALIZE VARIABLES
# =========================================================================

# Length of line segment between protein subdomains and exons
subdomain_gap <- 20

# Define colors to use for each subdomain and exon
protein_domain_colors <- c(
    "#AEC6CF", 
    "#CBAACB", 
    "#FFB7B2", 
    "#FF6961", 
    "#FFB347", 
    "#FDFD96", 
    "#77DD77", 
    "#B2FFFF"
)
exon_colors <- c(
    rep(protein_domain_colors[1], 15), # 15 exons corresponding to the first subdomain
    rep(protein_domain_colors[2], 3),  # 3 exons corresponding to the second subdomain
    rep(protein_domain_colors[3], 10), # 10 exons corresponding to the third subdomain
    rep(protein_domain_colors[4], 4),  # 4 exons corresponding to the fourth subdomain
    rep(protein_domain_colors[5], 2),  # 2 exons corresponding to the fifth subdomain
    rep(protein_domain_colors[6], 4),  # 4 exons corresponding to the sixth subdomain
    rep(protein_domain_colors[7], 6),  # 6 exons corresponding to the seventh subdomain
    rep(protein_domain_colors[8], 7)   # 7 exons corresponding to the eighth subdomain
)

# Start positions for each protein subdomain and exon
protein_domain_positions <- c(1, 705, 800, 1335, 1511, 1674, 1879, 2142, 2498)
exon_positions <- c(
    1, 287, 373, 483, 572, 707, 842, 974, 1094, 1237, 1317, 1424, 1554, 1679, 
    1792, 1937, 2077, 2206, 2377, 2636, 2825, 2944, 3014, 3232, 3483, 3632, 
    3726, 3913, 4095, 4325, 4453, 4672, 4874, 4963, 5151, 5306, 5453, 5645,
    5792, 5893, 6084, 6245, 6416, 6517, 6712, 6906, 6979, 7164, 7317, 7526, 
    7598, 9239
)


# =========================================================================
# LOAD GENOMIC VARIANT DATA TABLES
# =========================================================================

# Load data for each ancestry separately and also combined
df <- fread("lrrk2_data.tsv")
ancestry_tables <- split(df, df$Ancestry)
combined_table <- fread("lrrk2_grouped.tsv")
all_tables <- c(list(Combined = combined_table), ancestry_tables)

# Define function to process tables
clean_variant_table <- function(tbl) {
    # Find frequency of cases and controls
    tbl <- as.data.table(tbl)
    if (all(c("het_PD", "hom_PD", "total_PD") %in% names(tbl))) {
        tbl[, `PD Frequency` := (het_PD + hom_PD) / total_PD]
    }
    if (all(c("het_HC", "hom_HC", "total_HC") %in% names(tbl))) {
        tbl[, `Control Frequency` := (het_HC + hom_HC) / total_HC]
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
        "Mean_pRAB10/RAB10",
        "SD",
        "Interpretation"
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

    # Rename columns
    colnames(tbl) <- c(
        "Variant",
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
        "Kinase activity (mean pRAB10/RAB10)",
        "Standard deviation",
        "Interpretation"
    )

    return(tbl)
}

# Process tables, keeping combined tables by default
all_tables_cleaned <- lapply(all_tables, clean_variant_table)


# =========================================================================
# SERVER FUNCTION
# =========================================================================

server <- function(input, output, session) {
    # General metadata about LRRK2
    geneOverviewServer("gene_overview")

    # Counts of variants across functional annotation categories 
    annotationSummaryTableServer("annotation_summary_table", all_tables_cleaned)

    # Protein domain diagram
    proteinDiagramServer("protein_diagram", protein_domain_colors, protein_domain_positions, subdomain_gap)

    # cDNA diagram
    cdnaDiagramServer("cdna_diagram", exon_colors, exon_positions, subdomain_gap)

    # Main variant table
    geneVarTableServer("gene_var_table", all_tables_cleaned)
}
