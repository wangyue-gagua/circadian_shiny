singleCellUI <- function(id) {
  ns <- NS(id)
  
  tagList(plotOutput(ns("umap")),
          plotOutput(ns('violin')),
          plotOutput(ns('feature')))
}

singleCellServer <- function(id, plot_event, gene_id) {
  moduleServer(id,
               function(input, output, session) {
                 observeEvent(plot_event(), {
                   if (!str_detect(gene_id(),
                                   "^Ghir_\\w\\d{2}G\\d{6}$|^Ghir_\\w\\d{2}G\\d{6}\\.\\d$")) {
                     showNotification("Invalid input id!", type = "error")
                   } else{
                     ## transform from xx_xx to xx-xx term
                     id <- str_replace(gene_id(), '_', '-')
                     
                     output$umap <- renderPlot({
                       # generate bins based on input$bins from ui.R
                       DimPlot(
                         mergeSCE,
                         # reduction = "umap",
                         # pca, umap, tsne
                         group.by  = "Clu",
                         label = T
                       )
                     })
                     output$violin <- renderPlot({
                       VlnPlot(object = mergeSCE,
                               features = c(id))
                     })
                     output$feature <- renderPlot({
                       FeaturePlot(
                         mergeSCE,
                         # reduction = "umap",
                         features = c(id),
                         # split.by = "sample",
                         label = TRUE
                       )
                     })
                   }
                   
                   
                 })
               })
}