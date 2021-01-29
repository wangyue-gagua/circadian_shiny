

library(shiny)
library(tidyverse)
library(readxl)
library(shinythemes)
# library(gridExtra)


# import sample list , raw gene count and TMM normalized expression matrix
# genes.counts <- read.delim("./merged_counts/genes.counts.matrix", row.names=1, check.names = FALSE )
# genes.TMM.EXPR <- read.delim("./merged_counts/genes.TMM.EXPR.matrix", row.names=1, check.names = FALSE)


# data preprocess____________________________________________________________


# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("simplex"),
  # Application title
  titlePanel("Circadian rhythm plot generator"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      helpText(
        "Input a gene ID, e.g",
        em("Ghir_D11G029140"),
        br(),
        'click plot for expression, ',
        'GO information and homologs founded in TAIR10 '
      ),
      textInput(
        "geneid",
        label = h3("Input a gene ID"),
        placeholder = "Ghir_D11G029140"
      ),
      actionButton('plot', 'Plot'),
      downloadButton("downloadfig", "Download png"),
      width = 3
    ),
    
    mainPanel(tabsetPanel(
      tabPanel(
        "plots",
        h4("circadian rhythm, gene expression level"),
        plotOutput("circa_plot"),
        h4("expression level of gene in different tissues"),
        plotOutput("tissue_plot"),
        textOutput("tissue_err"),
        h4("different expression level of protein between FL and WT"),
        plotOutput("prot_plot"),
        textOutput("prot_err")
      ),
      tabPanel(
        "tables",
        h4("Go information"),
        DT::dataTableOutput('go'),
        h4("homologs in ", em("arabidopsis thaliana")),
        DT::dataTableOutput('homo_ara')
      )
    ))
    
    # Show a plot of the generated distribution
  )
)

#Define server logic required to draw a histogram
server <- function(input, output) {
  require(DT)
  
  observeEvent(input$plot,
               {
                 cir_plot <- my_cir_plot(input$geneid)
                 output$circa_plot <- renderPlot(cir_plot)
                 
                 tissue_plot <- my_tissue_explot(input$geneid)
                 if (is_character(tissue_plot)) {
                   output$tissue_err <- renderText(h1(tissue_plot))
                 } else{
                   output$tissue_plot <- renderPlot(tissue_plot)
                 }
                 
                 prot_plot <- my_prot_plot(input$geneid)
                 if (is_character(prot_plot)) {
                   output$prot_err <- renderText(h1(prot_plot))
                 } else{
                   output$prot_plot <- renderPlot(prot_plot)
                 }
                 
                 go_table <- my_go_table(input$geneid)
                 output$go <- DT::renderDataTable(go_table)
                 
                 df <- my_ara_homo_tbl(input$geneid)
                 if (!is.null(df$Match))
                   df <- df %>%
                   mutate(Match = map(Match,  ~ as.character(a(
                     href = str_c(
                       "https://www.arabidopsis.org/servlets/TairObject?type=locus&name=",
                       str_match(.x, "(^.*)\\.")[, 2]
                     ),
                     target = "_blank",
                     .x
                   ))))
                 output$homo_ara <-
                   DT::renderDataTable(expr = datatable(df, escape = FALSE))
               })
  
  # Downloadable csv of selected dataset ----
  output$downloadfig <- downloadHandler(
    filename = function() {
      paste(input$geneid, ".png", sep = "")
    },
    content = function(file) {
      ggsave(filename = file,
             plot = my_cir_plot(input$geneid),
             dpi = 500)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
