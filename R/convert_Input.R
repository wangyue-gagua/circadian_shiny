IdConvert_UI <- function(id) {
  ns <- NS(id)
  tagList(textInput(ns("ID"), h2("input an gene ID")),
          selectInput(
            ns("Genome"),
            h2("select your desired genome"),
            list(`G. hirsutum` = list("WHU", "ZJU"))
          ),
          actionButton(ns("Run_Convert"), "RUN"),)
}

IdConvert_Server <- function(id) {
  moduleServer(id, function(input, output, session) {
    return(list(
      id = reactive({
        input$ID
      }),
      genome = reactive({
        input$Genome
      }),
      run = reactive({
        input$Run_Convert
      })
    ))
  })
}