#!/usr/bin/env Rscript

# =========================================================================
# INITIALIZE VARIABLES
# =========================================================================

# Length of line segment between protein subdomains and exons
subdomain_gap <- 20

# Thresholds for kinase activity, CADD deleteriousness, and conservation
cadd_deleterious_threshold <- 20
conservation_conserved_threshold <- 5
kinase_activation_threshold <- 1.40
kinase_inactivation_threshold <- 0.50

# Define colors to use for each subdomain and exon
protein_domain_colors <- c(
    "#035C81",
    "#E3E3E3",
    "#0C8DC3",
    "#CCCCCC",
    "#4BA3C9",
    "#999999",
    "#86C1DA",
    "#666666"
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

# Store metadata for protein subdomains and exons
num_protein_domains <- length(protein_domain_positions) - 1
protein_domain_names <- c("ARM","ANK","LRR","ROC","COR-A","COR-B","KIN","WD-40")
protein_domain_positions_adj <- protein_domain_positions + (0:num_protein_domains) * subdomain_gap
protein_domain_tooltips <- c(
    "Armadillo repeats",
    "Ankyrin repeats",
    "Leucine-rich repeats",
    "Ras of complex protein",
    "C-terminal of ROC (A)",
    "C-terminal of ROC (B)",
    "Kinase",
    "WD-40"
)
protein_domains <- data.frame(
    start = protein_domain_positions_adj[-(num_protein_domains + 1)],
    end   = protein_domain_positions_adj[-1],
    label = protein_domain_names,
    color = protein_domain_colors,
    info  = protein_domain_tooltips
)
num_exons <- length(exon_positions) - 1
exon_names <- paste("", 1:num_exons)
exon_positions_adj <- exon_positions + (0:(num_exons)) * subdomain_gap
exon_tooltips <- paste("Exon", 1:num_exons)
exons <- data.frame(
    start = exon_positions_adj[-(num_exons+1)],
    end   = exon_positions_adj[-1],
    label = exon_names,
    color = exon_colors,
    info  = exon_tooltips
)


# =========================================================================
# LOAD GENOMIC VARIANT DATA TABLES
# =========================================================================

# Load data for each ancestry separately and also combined
df_imputed <- fread("lrrk2_combined_imputed.tsv")
df_wgs <- fread("lrrk2_combined_wgs.tsv")
df_raw <- fread("lrrk2_combined_raw.tsv")
df_exome <- fread("lrrk2_combined_exome.tsv")
all_tables_imputed <- split(df_imputed, df_imputed$Ancestry)
all_tables_wgs <- split(df_wgs, df_wgs$Ancestry)
all_tables_raw <- split(df_raw, df_raw$Ancestry)
all_tables_exome <- split(df_exome, df_exome$Ancestry)

# Process tables, keeping combined tables by default
all_tables_imputed_cleaned <- lapply(names(all_tables_imputed), function(name) {
    clean_variant_table(
        all_tables_imputed[[name]],
        ancestry = name,
        cadd_deleterious_threshold,
        conservation_conserved_threshold,
        kinase_activation_threshold,
        "Imputed"
    )
})
all_tables_wgs_cleaned <- lapply(names(all_tables_wgs), function(name) {
    clean_variant_table(
        all_tables_wgs[[name]],
        ancestry = name,
        cadd_deleterious_threshold,
        conservation_conserved_threshold,
        kinase_activation_threshold,
        "WGS"
    )
})
all_tables_raw_cleaned <- lapply(names(all_tables_raw), function(name) {
    clean_variant_table(
        all_tables_raw[[name]],
        ancestry = name,
        cadd_deleterious_threshold,
        conservation_conserved_threshold,
        kinase_activation_threshold,
        "Raw genotyping"
    )
})
all_tables_exome_cleaned <- lapply(names(all_tables_exome), function(name) {
    clean_variant_table(
        all_tables_exome[[name]],
        ancestry = name,
        cadd_deleterious_threshold,
        conservation_conserved_threshold,
        kinase_activation_threshold,
        "Clinical exome"
    )
})
names(all_tables_imputed_cleaned) <- names(all_tables_imputed)
names(all_tables_wgs_cleaned) <- names(all_tables_wgs)
names(all_tables_raw_cleaned) <- names(all_tables_raw)
names(all_tables_exome_cleaned) <- names(all_tables_exome)

# Merge imputed and WGS tables for each ancestry
all_tables_merged <- lapply(names(all_tables_imputed_cleaned), function(name) {
    imputed_table <- all_tables_imputed_cleaned[[name]]
    wgs_table <- all_tables_wgs_cleaned[[name]]
    raw_table <- all_tables_raw_cleaned[[name]]
    exome_table <- all_tables_exome_cleaned[[name]]
    
    # Perform outer merge on "Variant (GrCh38)"
    merged_table <- merge(
        imputed_table, 
        wgs_table, 
        by = "Variant (GrCh38)", 
        all = TRUE,
        suffixes = c("", ".wgs")
    )
    merged_table <- merge(
        merged_table, 
        raw_table,
        by = "Variant (GrCh38)", 
        all = TRUE, 
        suffixes = c("", ".raw")
    )
    merged_table <- merge(
        merged_table, 
        exome_table,
        by = "Variant (GrCh38)", 
        all = TRUE, 
        suffixes = c("", ".exome")
    )

    # For each base column name (strip .wgs/.raw/.exome), coalesce in priority order
    suffix_pat <- "(\\.wgs|\\.raw|\\.exome)$"
    base_cols <- unique(gsub(suffix_pat, "", setdiff(names(merged_table), "Variant (GrCh38)")))

    for (b in base_cols) {
        # columns available for this base, in priority order
        cand <- c(b, paste0(b, c(".wgs", ".raw", ".exome")))
        cand <- cand[cand %in% names(merged_table)]
        if (length(cand) >= 2) {
            merged_table[, (b) := do.call(fcoalesce, .SD), .SDcols = cand]
        }
    }

    # Drop the suffixed columns now that the base columns are filled
    drop_cols <- grep(suffix_pat, names(merged_table), value = TRUE)
    if (length(drop_cols)) merged_table[, (drop_cols) := NULL]

    merged_table <- unique(merged_table, by = "Variant (GrCh38)")
    
    return(merged_table)
})
names(all_tables_merged) <- names(all_tables_imputed_cleaned)

# Add a column indicating which variants are deemed pathogenic by GP2
pathogenic_variants <- fread("pathogenic_variants.txt", header = FALSE)$V1
all_tables_merged <- lapply(all_tables_merged, function(merged_table) {
    merged_table[, Pathogenic := ifelse(`Variant (GrCh38)` %in% pathogenic_variants, 1, 0)]
    return(merged_table)
})

# Variants to include in lollipops for domain diagrams
variants <- data.frame(
    aa_pos = c(1067, 1437, 1437, 1441, 1441, 1441, 1441, 1628, 1795, 2019, 2020, 2385),
    aa_label  = c("R1067Q","N1437H","N1437D","R1441S","R1441G","R1441C","R1441H","R1628P","L1795F","G2019S","I2020T","G2385R"),
    cdna_pos = c(3200, 4309, 4309, 4321, 4321, 4321, 4322, 4883, 5385, 6055, 6059, 7153),
    cdna_label  = c("C3200A","A4309C","A4309G","C4321A","C4321G","C4321T","G4322A","G4883C","G5385T","G6055A","T6059C","G7153A"),
    color  = c("#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444","#444444"),
    value  = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)
protein_variants <- variants[, c("aa_pos", "aa_label", "color", "value")]
cdna_variants <- variants[, c("cdna_pos", "cdna_label", "color", "value")]
colnames(protein_variants) <- c("pos", "label", "color", "value")
colnames(cdna_variants) <- c("pos", "label", "color", "value")


# =========================================================================
# SERVER FUNCTION
# =========================================================================

server <- function(input, output, session) {
    # Trigger disclaimer popup once when the session starts
    shinyalert(
        title = "Terms of Use",
        text = HTML("
            <div style='color: black; text-align: justify;'>
            <p style='color: black;'>By proceeding, I am agreeing to:</p>
            <ul>
                <li>Use the data for health, biomedical, and research ONLY</li>
                <li>NOT attempt to identify, disclose, or contact research participants 
                unless required by federal, state, or local laws</li>
                <li>Report any data management issues incidents, but not limited to, 
                inadvertent data release</li>
                <li>Abide by all relevant laws and regulations regarding genomic data 
                and their use</li>
                <li>NOT bulk download data without explicit consent from the LRRK2 
                Browser team </li>
            </ul>
            <p style='color: black;'>While data in this browser have undergone quality 
            control, genomic data processing pipelines are inherently imperfect and 
            rely on probabilistic processes such as variant calling, imputation, and 
            short-read sequencing reads. As such, there may be errors within the data 
            presented. I agree that the LRRK2 Browser team is not responsible for any 
            incorrect data that may be present in the browser.</p>
            </div>
        "),
        html = TRUE,
        type = "info"
    )
    
    # General metadata about LRRK2
    geneOverviewServer("gene_overview")

    # Links to other resources
    otherResourcesServer("other_resources")

    # Counts of variants across functional annotation categories 
    annotationSummaryTableServer("annotation_summary_table", all_tables_merged)

    # Protein domain diagram
    diagramServer("protein_diagram", protein_domains, protein_domain_positions, subdomain_gap, protein_variants, "protein", 12)

    # cDNA diagram
    diagramServer("cdna_diagram", exons, exon_positions, subdomain_gap, cdna_variants, "cDNA", 12)

    # Kinase activity bar chart
    exon_color_mapping <- setNames(
        rep(protein_domain_colors, times = c(17, 2, 9, 3, 3, 4, 5, 8)), 
        seq_len(length(exon_colors))
    )
    barChartServer("bar_chart", all_tables_merged$Combined, kinase_activation_threshold, kinase_inactivation_threshold, exon_color_mapping)

    # Main variant table with popup
    clicked_variant <- reactiveVal(NULL)
    geneVarTableServer("gene_var_table", all_tables_merged, clicked_variant)
    variantDetailServer("variant_detail", all_tables_merged, clicked_variant)
}
