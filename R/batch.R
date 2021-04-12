batchUI <- function(id, label = 'batch for heatmap') {
  ns <- NS(id)
  
  tagList(
    fileInput(ns('file_load'), 'gene id in tab-delimited format', accept = 'text/plain'),
    ### batch tab is active
    uiOutput(ns('pltHtMap')),
    h4("heatmap for multiple genes !"),
    strong('press F12 or resize your browser and then plots will appear'),
    plotOutput(ns("circa")),
    plotOutput(ns("tissue"))
  )
}

batchServer <- function(id) {
  moduleServer(
    id,
    
    function(input, output, session){
      observeEvent(input$file_load,{
        output$pltHtMap <- renderUI({
          ns <- session$ns
          actionButton(ns('heat_plot'), 'HeatMap')})
        
        id_list <-
          reactive({read_lines(input$file_load$datapath) %>% str_split('\t') %>% as_vector()})
        
        if_show_rownames = if_else(length(id_list()) < 20, true = TRUE, false = FALSE)
        
        df_circa <- genes.TMM.EXPR[id_list(), ]
        df_tissue <- iso_exp_tpm %>%
          filter(str_match(tracking_id, "(.*)\\.\\d")[,2] %in% id_list()) %>%
          column_to_rownames(var = "tracking_id")
        
        
          
          circa_plot <- reactive({
            pheatmap(
              df_circa %>% na.omit(),
              cluster_cols = FALSE,
              show_rownames = if_show_rownames,
              annotation_col = sample_info %>% column_to_rownames(var =
                                                                    'sample'),
              scale = 'none'
            )
          })
          
          output$circa <- renderPlot({
            circa_plot()
          })
          
          tissue_plot <- reactive({
            pheatmap(
              df_tissue %>% na.omit(),
              cluster_cols = FALSE,
              show_rownames = if_show_rownames,
              scale = 'none'
            )
          })
          
          output$tissue <- renderPlot({
            tissue_plot()
          })
        })
    }
  )
}