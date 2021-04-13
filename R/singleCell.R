singleCellUI <- function(id) {
  ns <- shiny::NS(id)

  tagList(
    div(
      class = "dim_redu_flex",
      img(src = "src/fiber_umap.png", width = "30vw"),
      img(src = "src/fiber_pca.png", width = "30vw"),
    ),
    div(
      class = "dim_redu_flex",
      img(src = "src/cor_gene_UMAP.png"),
      plotOutput(ns("sample_exp"))
    ),
    # plotOutput(ns("violin")),
    plotOutput(ns("feature")),
    plotOutput(ns("sample_vln"))
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
          # output$violin <- renderPlot({
          #   VlnPlot(
          #     object = mergeSCE,
          #     features = c(id)
          #   )
          # })
          output$sample_exp <- renderPlot({
            plotReducedDim(all_strain_sce,
              dimred = "UMAP",
              colour_by = gene_id()
            ) +
              facet_grid(. ~ all_strain_sce@colData$SampleName) +
              scale_color_gradient(low = "#cecbcb00", high = "#774aff")
          })
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