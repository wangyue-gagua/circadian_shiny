## for input id
geneInputUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    helpText(
      "Input a gene ID, e.g",
      em("Ghir_D11G029140, Ghir_D11G029140.1"),
      br(),
      'click plot for expression, ',
      'GO information and homologs found in TAIR10 '
    ),
    textInput(
      ns("geneid"),
      label = h3("Input a gene ID"),
      placeholder = "Ghir_D11G029140"
    ),
    #button
    actionButton(ns('plot'), 'Plot')
  )
}

geneInputServer <- function(id) {
  moduleServer(id,
               function(input, output, session) {
                 # 返回绘图事件和基因id
                 return(list(id = reactive({
                   input$geneid
                 }),
                 tri = reactive({
                   input$plot
                 })))
               })
}
