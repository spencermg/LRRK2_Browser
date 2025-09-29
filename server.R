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
    rep(protein_domain_colors[1], 17), # 17 exons corresponding to the first subdomain
    rep(protein_domain_colors[2], 2),  # 2 exons corresponding to the second subdomain
    rep(protein_domain_colors[3], 9),  # 9 exons corresponding to the third subdomain
    rep(protein_domain_colors[4], 3),  # 3 exons corresponding to the fourth subdomain
    rep(protein_domain_colors[5], 3),  # 3 exons corresponding to the fifth subdomain
    rep(protein_domain_colors[6], 4),  # 4 exons corresponding to the sixth subdomain
    rep(protein_domain_colors[7], 5),  # 5 exons corresponding to the seventh subdomain
    rep(protein_domain_colors[8], 8)   # 8 exons corresponding to the eighth subdomain
)

# Start positions for each protein subdomain and exon
protein_domain_positions <- c(1, 705, 800, 1335, 1511, 1674, 1879, 2142, 2498)
exon_positions <- c(
    1, 152, 238, 348, 437, 572, 707, 839, 959, 1102, 1182, 1289, 1419, 1544, 
    1657, 1802, 1942, 2071, 2242, 2501, 2690, 2809, 2879, 3097, 3348, 3497, 
    3591, 3778, 3960, 4190, 4318, 4537, 4739, 4828, 5016, 5171, 5318, 5510, 
    5657, 5758, 5949, 6110, 6281, 6382, 6577, 6771, 6844, 7029, 7182, 7391, 
    7463, 7585
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

# Variants to include in lollipops for the protein diagram
variants <- data.frame(
  aa_pos = c(1067, 1437, 1437, 1441, 1441, 1441, 1441, 1628, 1795, 2019, 2020, 2385),
  label  = c("R1067Q","N1437H","N1437D","R1441S","R1441G","R1441C","R1441H","R1628P","L1795F","G2019S","I2020T","G2385R"),
  color  = c("#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444"),
  value  = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

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
    proteinDiagramServer("protein_diagram", protein_domain_colors, protein_domain_positions, subdomain_gap, variants, 12)

    # cDNA diagram
    cdnaDiagramServer("cdna_diagram", exon_colors, exon_positions, subdomain_gap)

    # Main variant table
    geneVarTableServer("gene_var_table", all_tables_cleaned)
}
