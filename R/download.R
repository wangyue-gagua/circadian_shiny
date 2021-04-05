downloadUI <- function(id, label = 'download figure') {
  # namespace
  ns <- NS(id)
  
  tagList(
    downloadButton(ns("SourceData"), "Download data"),
    downloadButton(ns("fig"), "Download fig"),
    selectInput(
      ns("fileType"),
      label = "picture format",
      choices = list(
        "png" = ".png",
        "jpeg" = ".jpeg",
        "svg" = ".svg",
        "pdf" = ".pdf"
      )
    ),
  )
}


downloadServer <- function(id , gene_id) {
  moduleServer(
    id,
    ## module function
    function(input, output, session) {
      # Downloadable csv of selected dataset
      output$downloadfig <- downloadHandler(
        filename = function() {
          paste(gene_id, input$fileType, sep = "")
        },
        content = function(file) {
          ggsave(filename = file,
                 plot = my_cir_plot(gene_id))
        }
      )
      
      output$downloadSourceData <- downloadHandler(
        filename = function() {
          paste(gene_id, ".csv", sep = "")
        },
        content = function(filename) {
          sample_info_exp %>% select(sample, strain, period, time, replicate, labs, gene_id) %>%
            write_csv(file = filename, col_names = TRUE)
        }
      )
    }
  )
}