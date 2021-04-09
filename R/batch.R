batchUI <- function(id, label = 'batch for heatmap') {
  ns <- NS(id)
  
  tagList(
    fileInput(ns('file_load'), 'gene id in tab-delimited format', accept = 'text/plain'),
    ### batch tab is active
    uiOutput(ns('pltHtMap')),

  )
}

batchServer <- function(id) {
  moduleServer(
    id,
    
    function(input, output, session){
      observeEvent(input$file_load,{
        print('load')
        output$pltHtMap <- renderUI({
          ns <- session$ns
          actionButton(ns('heat_plot'), 'HeatMap')})
      })
      
      
      heatmap_RNA <- reactive({
        id_list <-
          read_lines(input$file_load$datapath) %>% str_split('\t') %>% as_vector()
        if_show_rownames = if_else(length(id_list) < 20, true = TRUE, false = FALSE)
        df <- genes.TMM.EXPR[id_list, ]
        pheatmap(
          df %>% na.omit(),
          cluster_cols = FALSE,
          show_rownames = if_show_rownames,
          annotation_col = sample_info %>% column_to_rownames(var =
                                                                'sample'),
          scale = 'none'
        )
      }) %>%
        bindCache(input$file_load$datapath) %>%
        bindEvent(input$heat_plot)
      
      return(heatmap_RNA)
    }
  )
}