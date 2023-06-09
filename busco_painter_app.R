if(!require(shiny)) install.packages("shiny", repos = "http://cran.us.r-project.org")
if(!require(shinydashboard)) install.packages("shinydashboard", repos = "http://cran.us.r-project.org")
if(!require(shinythemes)) install.packages("shinythemes", repos = "http://cran.us.r-project.org")
if(!require(ggplot2)) install.packages("shinythemes", repos = "http://cran.us.r-project.org")
if(!require(tidyverse)) install.packages("shinythemes", repos = "http://cran.us.r-project.org")
if(!require(scales)) install.packages("shinythemes", repos = "http://cran.us.r-project.org")
if(!require(rsconnect)) install.packages("shinythemes", repos = "http://cran.us.r-project.org")

#setwd('/Users/cw22/Documents/R_work/busco_painter_shiny_app/')

read_reference <- function(lineage){
    busco <- read_tsv(paste0(lineage, '_assignments.tsv'),
                      col_names = c("busco_id", "Status"),
                      col_types = c("cc"),
                      comment = "#") %>%
        return(busco)
}

read_busco <- function(buscoFile, reference_buscos, minimum){
    busco <- read.csv(buscoFile, sep='\t', skip=2, col.names =  c("busco_id", "Status", "Sequence",
                                                                                   "start", "end", "strand", "Score", "Length",
                                                                                   "OrthoDB_url", "Description")) %>%
        filter(Status == "Complete") %>%
        mutate(position = (start+end)/2) %>%
        select(busco_id, Sequence, position)
    busco <- merge(busco, reference_buscos, by="busco_id")[,c(1:4)]
    colnames(busco) <- c('buscoID', 'query_chr', 'position', 'assigned_element') 
    
    busco <- busco %>% group_by(query_chr) %>% 
        mutate(length = max(position)) %>% 
        mutate(n_busco = n()) %>% # make a new column reporting number buscos per query_chr
        ungroup() %>%
        mutate(start = 0) %>%
        filter(n_busco >= minimum)  # filter df to only keep query_chr with >=3 buscos to remove shrapnel
    
    return(busco)
}


busco_paint_theme <- theme(legend.position="right",
                           strip.text.x = element_text(margin = margin(0,0,0,0, "cm")), 
                           panel.background = element_rect(fill = "white", colour = "white"), 
                           panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(),
                           axis.line.x = element_line(color="black", size = 0.5),
                           axis.text.x = element_text(size=10),
                           axis.title.x = element_text(size=10),
                           strip.text.y = element_text(angle=0),
                           strip.background = element_blank(),
                           plot.title = element_text(hjust = 0.5, face="italic", size=15),
                           plot.subtitle = element_text(hjust = 0.5, size=12))

busco_paint <- function(spp_df, num_col, title, karyotype, ALG_name){
  #  colour_palette <- append(hue_pal()(32), 'grey')
   # merian_order <- c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
    colour_palette <- hue_pal()(length(unique(spp_df$assigned_element)))
   # spp_df$assigned_element_f =factor(spp_df$assigned_element, levels=merian_order)
    spp_df$assigned_element_f =factor(spp_df$assigned_element, levels=unique(spp_df$assigned_element))
    chr_levels <- subset(spp_df, select = c(query_chr, length)) %>% unique() %>% arrange(length, decreasing=TRUE)
    chr_levels <- chr_levels$query_chr
    spp_df$query_chr_f =factor(spp_df$query_chr, levels=chr_levels) # set chr order as order for plotting
    sub_title <- paste("n contigs =", karyotype) 
    the_plot <- ggplot(data = spp_df) +
        scale_colour_manual(values=colour_palette, aesthetics=c("colour", "fill")) +
        geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
        geom_rect(aes(xmin=position-2e4, xmax=position+2e4, ymax=0, ymin =12, fill=assigned_element_f)) +
        facet_wrap(query_chr_f ~., ncol=num_col, strip.position="right") + guides(scale="none") + 
        xlab("Position (Mb)") +
        scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
        scale_y_continuous(breaks=NULL) + 
        ggtitle(label=title, subtitle= sub_title)  +
        guides(fill=guide_legend(ALG_name), color = "none") +
        busco_paint_theme
    return(the_plot)
}

# inputs before running server


# busco painter, about this site
### specify UI

ui <-
    navbarPage(theme = shinytheme("flatly"), collapsible = TRUE,
               HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">BUSCO painter</a>'), id="nav",
               windowTitle = "BUSCO painter",
               tabPanel("BUSCO paint",
                        fluidPage(
                                sidebarLayout(
   
                                    sidebarPanel(
                                        selectInput("Lineage", "Select the orthoDB dataset used", c("Lepidoptera_odb10", "Nematoda_odb10")),
                                        selectInput("ALG", "Select the ancestral linkage groups of interest", c("Merian elements", "Nigon units")), 
                                        fileInput("file1", "Upload a full busco table in tsv format", accept = c(".tsv")), # only accept tsv files
                                        numericInput("minimum", "Filter scaffolds with less than a minimum number of BUSCOs", value = 3, min = 1, max=20),
                                        numericInput("c", "Number of columns in plot", value = 2, min = 1, max=5),
                                        textInput("title", "Title of plot", placeholder="Title..", value=""), # placeholder=label inside the text box
                                    ),
                                    
                                    mainPanel(
                                        tabsetPanel(
                                            tabPanel("Plot", plotOutput("plot")), 
                                            tabPanel("Summary", tableOutput("summary")), 
                                            tabPanel("Table", tableOutput("table")),
                                        )
                                    )
                                ),
                                fluidRow(
                                    downloadButton("downloadPlot", "Download plot"),
                                    downloadButton("downloadTable", "Download tsv"))
               )
               ),
               tabPanel("About this site",
                                 tags$div(
                                   tags$h4("Last update"), 
                                   h6("7th June 2023"),
                                     "This site will be updated as new ancestral linkage group declarations become available.",
                                     tags$br(),tags$br(),tags$h4("Background"), 
                                     "Ancestral linkage groups (ALGs) are the groups of linked genes present in the last common ancestor of a given taxa. 
                                    The increasing availability of chromosomal genomes across taxa, mean that ALGs for groups of species can be inferred. 
                                   ALGs have been inferred for Lepidoptera (butterflies and moths), Diptera, flowering plants, Nematoda, mammals, vertebrates and Metazoa. 
                                   This tool allows a given genome to be compared to a chosen set of ALG by painting the positions of orthologs that define each ALG onto the chromosomes of the genome. This enables the evolutionary history of ALGs to be explored in the given genome. All that is required as input is the full table generated by running",
                                   tags$a(href="https://academic.oup.com/bioinformatics/article/31/19/3210/211866?login=false", "BUSCOs"), "on the genome of interest.",
                                     tags$br(),tags$br(),tags$h4("Code"),
                                   "Code and reference data used to generate this Shiny painting tool are available at ",tags$a(href="https://github.com/charlottewright/busco_painter", "Github."),
                                     tags$br(),tags$br(),tags$h4("Cite"),
                                  tags$b("Merian elements (Lepidoptera ALGs): "), tags$a(href="https://www.biorxiv.org/content/10.1101/2023.05.12.540473v1", "Chromosome evolution in Lepidoptera"), "CJ Wright, L Stevens, A Mackintosh, M Lawniczak, M Blaxter, bioRxiv, 2023.05 12.540473",tags$br(),
                                  tags$b("Nigon units (Rhabditid nematoda ALGs): "), tags$a(href="https://academic.oup.com/g3journal/article/11/1/jkaa020/6026964", "A telomere-to-telomere assembly of Oscheius tipulae and the evolution of rhabditid nematode chromosomes"), "P Gonzalez de la Rosa, M Thomson, U Trivedi, A Tracey, S Tandonnet, M Blaxter. G3, 11 (1), 2021",
                                     tags$br(),tags$br(),tags$h4("Contact"),
                                     "charlotte.wright@sanger.ac.uk",tags$br(),tags$br()
                                    # tags$img(src = "vac_dark.png", width = "150px", height = "75px"), tags$img(src = "lshtm_dark.png", width = "150px", height = "75px")
                                 )
                        )
    )         
    


## Specify backend                        
server <- function(input, output, session) {
    lineage2alg_name <- c("Lepidoptera_odb10"="Merian element", "Nematoda_odb10"="Nigon unit")
 #   ALG_name <- reactive({lineage2alg_name[input$Lineage]})
    assignments <- reactive({
        req(input$Lineage)
        assignments <- read_reference(input$Lineage)
        })
        data <- reactive({
            req(input$file1) # make sure the code waits until the first file is uploaded
            data <- read_busco(input$file1$datapath, assignments(), input$minimum)
    })

    num_contigs <- reactive({
        as.character(length(unique(data()$query_chr))) # number of query_chr after filtering
    })
    
    plotInput<- reactive({
        busco_paint(data(), input$c, input$title, num_contigs(), input$ALG)
    })
    
    output$plot <- renderPlot({ 
        plotInput()
        }, res=96)
    
    output$table <- renderTable({
        head(data()[,c(1:4)], 20)
    })
    
    output$summary <- renderTable({
        data() %>% group_by(query_chr, assigned_element) %>% mutate(n_markers = n()) %>% select(query_chr, assigned_element, n_markers) %>% arrange(desc(n_markers)) %>% distinct()
    })
    # Download plot
    output$downloadPlot <- downloadHandler(
        filename = function() { paste("TEST", '.png', sep='') },
        content = function(file) {
            ggsave(file, plot = plotInput(), device = "png")
        }
    )
    
    output$downloadTable <- downloadHandler(
        filename = function() {
            paste0(input$Lineage, ".paint", ".tsv")
        },
        content = function(file) {
            write.table(data()[,c(1:4)], file, sep='\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
        }
    )
}

### Run server
shinyApp(ui, server)

