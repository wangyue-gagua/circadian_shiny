singleCellUI <- function(id) {
  ns <- shiny::NS(id)

  tagList(
    div(
      class = "dim_redu_flex",
      img(src = "src/all_genes_UMAP.png", width = "30vw"),
      img(src = "src/cor_gene_UMAP.png", width = "30vw")
    ),
    div(
      class = "dim_redu_flex",
      plotOutput(ns("umap")),
      plotOutput(ns("pca"))
    ),
    # plotOutput(ns("violin")),
    plotOutput(ns("sample_vln")),
    plotOutput(ns("feature"))
  )
}

singleCellServer <- function(id, plot_event, gene_id) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(plot_event(), {
        if (!str_detect(
          gene_id(),
          "^Ghir_\\w\\d{2}G\\d{6}$|^Ghir_\\w\\d{2}G\\d{6}\\.\\d$"
        )) {
          shiny::showNotification("Invalid input id!", type = "error")
        } else {
          ## transform from xx_xx to xx-xx term
          id <- stringr::str_replace(gene_id(), "_", "-")

          output$umap <- renderPlot({
            # generate bins based on input$bins from ui.R
            DimPlot(
              mergeSCE,
              reduction = "UMAP",
              # pca, umap, tsne
              group.by  = "Clu",
              label = T
            )
          })
          output$pca <- renderPlot({
            DimPlot(mergeSCE,
              reduction = "PCA",
              group.by = "RawPCAClu", label = T
            )
          })
          # output$violin <- renderPlot({
          #   VlnPlot(
          #     object = mergeSCE,
          #     features = c(id)
          #   )
          # })
          output$sample_vln <- renderPlot({
            plotExpression(all_strain_sce, c(gene_id()),
              x = "cluster", colour_by = "SampleName"
            )
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
    }
  )
}