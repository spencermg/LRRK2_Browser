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
    "#035c81",
    "#e3e4e3",
    "#0c8dc3", 
    "#cccccc", 
    "#4ba3c9", 
    "#999999", 
    "#86c1da", 
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
protein_domain_tooltips <- rep("Protein subdomain", num_protein_domains)
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
df <- fread("lrrk2_data.tsv")
ancestry_tables <- split(df, df$Ancestry)
combined_table <- fread("lrrk2_grouped.tsv")
all_tables <- c(list(Combined = combined_table), ancestry_tables)

# Process tables, keeping combined tables by default
all_tables_cleaned <- lapply(all_tables, function(x) clean_variant_table(
    x, 
    cadd_deleterious_threshold,
    conservation_conserved_threshold,
    kinase_activation_threshold
))

### TODO: Update this with all variants and real values for variants
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
    # General metadata about LRRK2
    geneOverviewServer("gene_overview")

    # Counts of variants across functional annotation categories 
    annotationSummaryTableServer("annotation_summary_table", all_tables_cleaned)

    # Protein domain diagram
    diagramServer("protein_diagram", protein_domains, protein_domain_positions, subdomain_gap, protein_variants, "protein", 12)

    # cDNA diagram
    diagramServer("cdna_diagram", exons, exon_positions, subdomain_gap, cdna_variants, "cDNA", 12)

    # Kinase activity bar chart
    exon_color_mapping <- setNames(
        rep(protein_domain_colors, times = c(17, 2, 9, 3, 3, 4, 5, 8)), 
        seq_len(length(exon_colors))
    )
    barChartServer("bar_chart", all_tables_cleaned$Combined, kinase_activation_threshold, kinase_inactivation_threshold, exon_color_mapping)

    # Reactive value to store clicked variant information
    clicked_variant <- reactiveVal(NULL)

    # Main variant table
    geneVarTableServer("gene_var_table", all_tables_cleaned, clicked_variant)
    
    # Variant detail modal
    variantDetailServer("variant_detail", clicked_variant, all_tables_cleaned)
}
