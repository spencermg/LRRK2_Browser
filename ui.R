#!/usr/bin/env Rscript

# =========================================================================
# INITIALIZE VARIABLES
# =========================================================================

metadata_box_color         <- "#0C8DC3" # GP2 dark blue
variant_function_box_color <- "#D8EAF8" # GP2 light blue
protein_diagram_box_color  <- "#EAF4FF" # GP2 light grey
cdna_diagram_box_color     <- "#8C4E9F" # GP2 purple
variant_table_color        <- "#CF2FB3" # GP2 pink
box_outline_color          <- "#D3D3D3" # Light grey
box_title_font_style       <- "font-weight: bold; font-size: 28px;"


# =========================================================================
# UI FUNCTION
# =========================================================================

ui <- dashboardPage(
    title = "LRRK2 Browser",
    header = dashboardHeader(),
    sidebar = dashboardSidebar(
        width = 350,
        sidebarMenu(
            id = "sidebar"
        ),
        tags$div(
            style = "padding: 20px; overflow-y: auto; height: calc(100vh - 100px);",
            tags$h3("About the LRRK2 Browser", style = "color: #0C8DC3; margin-top: 0;"),
            tags$hr(),
            
            tags$h4("Data Sources"),
            tags$p(
                "This browser contains genomic data for participants from the",
                tags$a(
                    href = "https://doi.org/10.5281/zenodo.17753486", 
                    target = "_blank",
                    style = "color: #0C8DC3 !important;",
                    "11th data release"
                ),
                "of the Global Parkinson's Genetics Program (GP2). These data include several modalities:"
            ),
            tags$ul(
                tags$li(tags$strong("Whole genome sequencing (WGS)")),
                tags$li(tags$strong("Imputed genotyping")),
                tags$li(tags$strong("Raw genotyping")),
                tags$li(tags$strong("Clinical exome sequencing (CES)"))
            ),
            
            tags$h4("Sample Information"),
            tags$p("Total samples: 101,678"),
            tags$ul(
                tags$li("WGS: 25,904"),
                tags$li("Imputed genotyping: 84,654"),
                tags$li("Raw genotyping: 84,654"),
                tags$li("CES: 14,555")
            ),
            tags$p("Parkinson's Disease cases: 61,709"),
            tags$ul(
                tags$li("WGS: 17,241"),
                tags$li("Imputed genotyping: 46,349"),
                tags$li("Raw genotyping: 46,349"),
                tags$li("CES: 14,555")
            ),
            tags$p("Healthy controls: 39,969"),
            tags$ul(
                tags$li("WGS: 8,663"),
                tags$li("Imputed genotyping: 38,305"),
                tags$li("Raw genotyping: 38,305"),
                tags$li("CES: 0")
            ),
            
            tags$h4("Citation"),
            tags$p("For more information, or if you use data from this browser, please refer to:"),
            tags$p(
                tags$em("Citation"),
                style = "background-color: #f5f5f5; padding: 10px; border-left: 3px solid #0C8DC3;"
            ),
            
            tags$h4("Contact"),
            tags$p("For questions, please contact:"),
            tags$p(
                tags$a(
                    href = "mailto:lara.lange@nih.gov", "lara.lange@nih.gov",
                    style = "color: #0C8DC3 !important;"
                )
            ),
            tags$p(
                tags$a(
                    href = "mailto:spencer.grant@nih.gov", "spencer.grant@nih.gov",
                    style = "color: #0C8DC3 !important;"
                )
            ),
            tags$p(
                tags$a(
                    href = "mailto:vesna.van.midden@gmail.com", "vesna.van.midden@gmail.com",
                    style = "color: #0C8DC3 !important;"
                )
            ),
            

            tags$h3("About GP2", style = "color: #0C8DC3; margin-top: 50px;"),
            tags$hr(),

            tags$h4("General info"),
            tags$p(
                "The ",
                tags$a(
                    href = "https://gp2.org/", 
                    target = "_blank",
                    style = "color: #0C8DC3 !important;",
                    "Global Parkinson's Genetics Program"
                ),
                "(GP2) is a resource program of the",
                tags$a(
                    href = "https://parkinsonsroadmap.org/", 
                    target = "_blank",
                    style = "color: #0C8DC3 !important;",
                    "Aligning Science Across Parkinson's"
                ),
                "(ASAP) initiative, which is managed by the",
                tags$a(
                    href = "https://www.aligningscience.com/", 
                    target = "_blank",
                    style = "color: #0C8DC3 !important;",
                    "Coalition for Aligning Science"
                ),
                "(CAS) and implemented by the ",
                tags$a(
                    href = "https://www.michaeljfox.org/", 
                    target = "_blank",
                    style = "color: #0C8DC3 !important;",
                    "Michael J. Fox Foundation"
                ),
                paste0(
                    "(MJFF). GP2 serves as a central hub to generate genetic data from studies across the globe and ",
                    "harmonize corresponding clinical data. Our goal is to genetically characterize over 250,000 volunteers ",
                    "around the world to further understand the genetic architecture of Parkinson's disease (PD). For data ",
                    "processing and QC pipelines, please visit our "
                ),
                tags$a(
                    href = "https://github.com/GP2code/", 
                    target = "_blank",
                    style = "color: #0C8DC3 !important;",
                    "Github repository."
                )
            ),

            tags$h4("GP2 member sites and investigators"),
            tags$p(
                "A full list of GP2 sites and principal investigators can be found",
                tags$a(
                    href = "https://zenodo.org/records/17753486", 
                    target = "_blank",
                    style = "color: #0C8DC3 !important;",
                    "here."
                )
            ),

            tags$h4("Acknowledgments"),
            tags$p(paste0(
                "We would like to thank the many tens of thousands of participants who generously contributed to our effort. ",
                "In addition, we want to sincerely thank all the GP2 members who are involved in the program."
            )),

            
            tags$h3("Terms of Use", style = "color: #0C8DC3; margin-top: 50px;"),
            tags$hr(),

            tags$p(
                paste0(
                    "This website is intended to provide summary-level results from high-throughput sequencing data generated as ",
                    "part of GP2. Data and results are intended solely for educational and research purposes. The results, browser ",
                    "content and allele frequencies are NOT for diagnostic, clinical, or commercial use and must not be used to attempt ",
                    "to identify individual study participants.",
                    "We encourage you to contact GP2 before embarking on any analyses using these data to check if your proposed analyses ",
                    "overlap with work currently underway by GP2 members."
                )
            )
        )
    ),
    body = dashboardBody(
        # Initialize shinyjs
        shinyjs::useShinyjs(),
        
        # CSS to fix sidebar colors and hide content initially
        tags$head(
            tags$style(HTML("
                /* Disable bounce/overscroll effect */
                html, body {
                    overscroll-behavior-y: none !important;
                }
                
                /* Hide content until authenticated */
                .content-wrapper {
                    display: none;
                }

                .main-sidebar {
                    display: none;
                }
                
                /* Show when authenticated class is added */
                .content-wrapper.authenticated {
                    display: block !important;
                }
                
                .main-sidebar.authenticated {
                    display: block !important;
                }
                
                /* Blur/white out background when modal showing */
                .modal-backdrop {
                    backdrop-filter: blur(10px);
                    background-color: rgba(255, 255, 255, 0.95) !important;
                }
                
                /* Make entire header uniformly blue - ALL parts */
                .main-header {
                    background-color: #367fa9 !important;
                }
                
                .main-header .logo {
                    background-color: #367fa9 !important;
                }
                
                .main-header .navbar {
                    background-color: #367fa9 !important;
                    margin-left: 350px !important;
                }
                
                /* When sidebar is collapsed, no margin */
                .sidebar-collapse .main-header .navbar {
                    margin-left: 0px !important;
                }

                /* Style the toggle button */
                .sidebar-toggle {
                    color: white !important;
                    padding-right: 15px !important;
                }
                
                /* Add 'About the Data' label after the three lines */
                .sidebar-toggle::after {
                    content: 'About the Data';
                    margin-left: 10px;
                    font-size: 14px;
                    font-weight: bold;
                    font-family: 'Arial', sans-serif;
                    font-style: normal;
                    letter-spacing: 0.5px;
                    text-transform: none;
                    vertical-align: middle;
                }
                
                /* Make entire sidebar background match header blue */
                .main-sidebar {
                    background-color: #367fa9 !important;
                    position: fixed !important;
                    top: 0 !important;
                    left: 0 !important;
                    height: 100vh !important;
                    overflow-y: auto !important;
                }
                
                /* The content area inside sidebar should be white/light */
                .main-sidebar .sidebar {
                    background-color: #ffffff !important;
                }
                
                /* Make all text dark/visible */
                .sidebar h3, .sidebar h4, .sidebar p, .sidebar li, .sidebar a {
                    color: #333333 !important;
                }
                
                /* Hide text content when sidebar is collapsed */
                .sidebar-collapse .main-sidebar .sidebar {
                    display: none !important;
                }
                
                /* Fix the header portion of sidebar to be blue */
                .main-sidebar .sidebar-menu {
                    background-color: #367fa9 !important;
                }
            ")),
    
            tags$script(HTML("
                Shiny.addCustomMessageHandler('show_content', function(message) {
                    $('.content-wrapper').addClass('authenticated');
                    $('.content-wrapper').css('display', 'block');
                    $('.main-sidebar').addClass('authenticated');
                    $('.main-sidebar').css('display', 'block');
                });
            "))
        ),
        
        # Title banner
        h2("LRRK2 Browser", style = "text-align:center; font-weight:bold;"),

        # Make collapse functionality more obvious
        tags$head(
            tags$style(HTML("
                /* Add text before the icon, and some spacing */
                .box .btn-box-tool[data-widget='collapse']::before {
                content: 'Collapse ';
                font-weight: bold;
                }

                /* When the box is collapsed, change the text to 'Expand' */
                .box.collapsed-box .btn-box-tool[data-widget='collapse']::before {
                content: 'Expand ';
                }
                
                /* Optional: Adjust padding if the text feels too cramped */
                .box .btn-box-tool {
                padding: 5px 10px;
                }
            "))
        ),

        # Overview numbers for LRRK2
        box(
            title = tags$div(
                "LRRK2 Metadata",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", metadata_box_color, "; padding: 10px;"),
            tags$div(
                "Leucine Rich Repeat Kinase 2",
                style = "text-align: center; font-weight: bold; font-size: 20px; margin-bottom: 10px;"
            ),
            geneOverviewUI("gene_overview")
        ),

        # Main variant table
        box(
            title = tags$div(
                "Variant Table",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", variant_table_color, "; padding: 10px;"),
            geneVarTableUI("gene_var_table")
        ),

        # Summary table for different variant functional annotation categories
        box(
            title = tags$div(
                "Variant Function Summary",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", variant_function_box_color, "; padding: 10px;"),
            annotationSummaryTableUI("annotation_summary_table")
        ),

        # Protein subdomain diagram
        box(
            title = tags$div(
                "Protein Diagram",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", protein_diagram_box_color, "; padding: 10px;"),
            diagramUI("protein_diagram")
        ),

        # cDNA diagram
        box(
            title = tags$div(
                "cDNA Diagram",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", cdna_diagram_box_color, "; padding: 10px;"),
            diagramUI("cdna_diagram")
        ),

        # cDNA diagram
        box(
            title = tags$div(
                "Kinase Activity by Variant",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", cdna_diagram_box_color, "; padding: 10px;"),
            barChartUI("bar_chart")
        ),

        # Links to external resources
        box(
            title = tags$div(
                "Other Resources",
                style = box_title_font_style
            ),
            closable = FALSE,
            collapsible = TRUE,
            collapsed = FALSE,
            status = NULL,
            width = 12,
            solidHeader = FALSE,
            style = paste0("border: 1px solid ", box_outline_color, "; border-top: 3px solid ", metadata_box_color, "; padding: 10px;"),
            otherResourcesUI("other_resources")
        )
    ),
    controlbar = NULL
)
