### GENE SET ENRICHMENT EXPLORER

libDir <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/4.1"
if (dir.exists(libDir)) .libPaths(libDir)

library(shiny)
library(tidyverse)
library(DT)
library(kableExtra)
source("../shiny_func.R")

select <- dplyr::select

css <- "table{font-size: 11px; background-color: #EAF5FF}"

### Read data ###

data <- sh_read_gsea_data("../data")
initial_contrasts <- unique(data$fg$sel$go$contrast) %>% as.character()


#######################################################################

ui <- shinyUI(fluidPage(
  
  #tags$style(css),
  
  titlePanel("Tfe1/Tfe2 toxin expression in yeast: gene set enrichment"),

  radioButtons("terms", "Functional terms", choices = c("GO" = "go", "Reactome" = "re", "KEGG" = "kg"), inline = TRUE),
  radioButtons("method", "Method:", choices = c("Pairwise" = "sel", "Full model" = "fi", "Tfe correlation" = "tfe"), inline = TRUE),
  selectInput("contrast", "Contrast:", choices = initial_contrasts),

  fluidRow(
    column(6, 
      dataTableOutput("term_table"),
      fluidRow(
        column(6, 
               br(),
               radioButtons("intensity_scale", "Intesity scale:", choices = c("lin" = "", "log" = "log"), inline = TRUE),
               plotOutput("gene_plot", height = "400px",width = "50%")
        ),
        column(6, br(), tableOutput("gene_info"))
      ),
    ),
    column(6, dataTableOutput("gene_table"))
  ),
  
  
  
))


########################################################################################


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Prevents RStudio from crashing when Shiny window closed manually
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # Available contrasts
  get_contrasts <- function() {
    data$fg[[input$method]][[input$terms]]$contrast %>% 
      as.character() %>% 
      unique()
  }
  
  # Update contrasts based on method
  observeEvent(input$method, {
    freezeReactiveValue(input, "contrast")
    updateSelectInput(session = session, inputId = "contrast", choices = get_contrasts())
  })
  
  
  get_data <- function() {
    ctr <- input$contrast
    data$fg[[input$method]][[input$terms]] %>% 
      filter(contrast == ctr)
  }
  
  get_term_data <- function() {
    fg <- get_data()
    fg %>% 
      filter(padj < 0.05)
  }
  
  get_de <- function(filter_contrast = TRUE) {
    mthd <- input$method
    ctrs <- input$contrast
    # in case of tfe use factorial de instead
    if (mthd == "tfe") {
      mthd <- "fi"
    }
    if (ctrs == "corTfe1") {
      ctrs <- "strainTfe1"
    } else if (ctrs == "corTfe2") {
      ctrs <- "strainTfe2"
    }
    d <- data$de[[mthd]]
    if (filter_contrast) d <- d %>% filter(contrast == ctrs)
    d
  }
  
  get_gene_data <- function(trm) {
    fg <- get_data()
    de <- get_de()
    print(de)
    r <- fg %>% 
      filter(term == trm)
    df <- NULL
    if (nrow(r) > 0) {
      genes <- r$leading_edge[[1]]
      df <- de %>% 
        filter(gene_id %in% genes)
    }
    df
  }
  
  select_term <- function() {
    sel <- NULL
    term_data <- get_term_data()
    tab_idx <- as.numeric(input$term_table_rows_selected)
    if (length(tab_idx) == 1) {
      sel <- term_data[tab_idx, ] %>% pull(term)
    }
    return(sel)
  }
  
  select_gene <- function() {
    sel <- NULL
    term <- select_term()
    if (!is.null(term)) {
      gene_data <- get_gene_data(term)
      tab_idx <- as.numeric(input$gene_table_rows_selected)
      if (length(tab_idx) > 0) {
        gd <- gene_data[tab_idx, ]
        if (!is.null(gd)) sel <- gd %>% pull(gene_id)
      }
    }
    return(sel)
  }

  output$term_table <- DT::renderDataTable({
    terms <- get_term_data() %>% 
      select(term, term_name, NES) %>% 
      mutate_if(is.numeric, ~signif(.x, 2)) %>% 
      mutate(term = as.character(term))
    terms
  }, selection = "single", escape = FALSE)
  
  
  output$gene_table <- DT::renderDataTable({
    term <- select_term()
    df <- NULL
    if (!is.null(term)) {
      gene_data <- get_gene_data(term)
      if (!is.null(gene_data)) {
        df <- gene_data %>% 
          select(gene_name, description, gene_biotype, logFC) %>% 
          mutate_if(is.numeric, ~signif(.x, 2))
      }
    }
    df
  }, selection = "single")


  output$gene_plot <- renderPlot({
    dat <- get_data()
    sel <- select_gene()
    if (!is.null(sel)) {
      sh_plot_genes(data$star, sel, input$intensity_scale)
    }
  })
  
  output$gene_info <- function() {
    de <- get_de(filter_contrast = FALSE)
    sel <- select_gene()
    if (!is.null(sel)) {
      d <- de %>% 
        filter(gene_id %in% sel) %>% 
        mutate_if(is.numeric, ~signif(.x, 2)) %>% 
        select(logFC, logCPM, PValue, FDR, contrast)
      d %>% 
        kable("html") %>% 
        kable_styling("striped", full_width = F) %>%
        row_spec(which(d$FDR < 0.05 & abs(d$logFC) > 1), bold = TRUE, background = "lightsalmon")
        
    }
  }
  
}

# Run the application
shinyApp(ui = ui, server = server)
