ConvertOut_UI <- function(id) {
  ns <- NS(id)
  tagList(tableOutput(ns("ID_convert")))
}

ConvertOutServer <- function(id, genomeVer, geneID, plotEvent) {
  moduleServer(id, function(input, output, session) {
    observeEvent(plotEvent(), {
      if (genomeVer() == "WHU") {
        # Input validation
        if (str_detect(geneID(), "[AD]\\d{2}G\\d+")) {
          result <- IdConvert(WHU, geneID())
          if (nrow(result) == 0)
            result <- "No data available"
        } else {
          result <- "Wrong ID Format"
        }
      }
      output$ID_convert <- renderTable({
        result
      })
    })
  })
}