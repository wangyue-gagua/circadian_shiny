



# data preprocess____________________________________________________________

### 封装失败,toggle功能正常,图片能正常生成,问题很有可能在于output,渲染UI的位置不对,但不知如何解决
# textPlotOutput <- function(fn_plot, id_plot, id_err) {
#   textPlot <- fn_plot
#   if (is_character(textPlot)) {
#     nsid <- str_c('plots',id_err,sep = '-')
#     output$nsid <- renderText({textPlot})
#     print(textPlot)
#   }else{
#     nsid <- str_c('plots',id_plot,sep = '-')
#     output$nsid <- renderPlot(textPlot)
#     print(str(textPlot))
#   }
#
#   toggleElement(id = id_plot,
#                 condition = !is_character(textPlot))
#   toggleElement(id = id_err, condition = is_character(textPlot))
# }

# data preprocess____________________________________________________________




# _____________________________________________________________________________

plotUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("circadian rhythm, gene expression level"),
    plotOutput(ns("circa_plot")),
    textOutput(ns("circa_err"), container = h1),
    h4("expression level of isoform in different tissues"),
    plotOutput(ns("tissue_plot")),
    textOutput(ns("tissue_err"), container = h1),
    h4("different expression level of protein between FL and WT"),
    plotOutput(ns("prot_plot")),
    textOutput(ns("prot_err"), container = h1)
  )
}

plotServer <- function(id, plot_event, gene_id) {
  moduleServer(
    id,
    function(input, output, session) {


      ###
      # default plot page
      de_cir <- my_cir_plot("Ghir_D11G029140")

      de_tissue <- my_tissue_plot("Ghir_D11G029140")
      de_pro <- my_prot_plot("Ghir_D11G029140")


      output$circa_plot <- renderPlot(de_cir)
      output$tissue_plot <-
        renderPlot(de_tissue)

      output$"prot_plot" <- renderPlot(de_pro)


      output$prot_plot <- renderPlot(de_pro)


      observeEvent(plot_event(), {
        if (!str_detect(
          gene_id(),
          "^Ghir_\\w\\d{2}G\\d{6}$|^Ghir_\\w\\d{2}G\\d{6}\\.\\d$"
        )) {
          showNotification("Invalid input id!", type = "error")
        } else {
          circa_plot <-
            my_cir_plot(str_extract(gene_id(), "^Ghir_\\w\\d{2}G\\d{6}"))

          if (is_character(circa_plot)) {
            output$circa_err <- renderText(circa_plot)
          } else {
            output$circa_plot <- renderPlot(circa_plot)
          }
          toggleElement(
            id = "circa_plot",
            condition = !is_character(circa_plot)
          )
          toggleElement(id = "circa_err", condition = is_character(circa_plot))

          tissue_plot <- my_tissue_plot(gene_id())
          if (is_character(tissue_plot)) {
            output$tissue_err <- renderText(tissue_plot)
          } else {
            output$tissue_plot <- renderPlot(tissue_plot)
          }
          toggleElement(
            id = "tissue_plot",
            condition = !is_character(tissue_plot)
          )
          toggleElement(id = "tissue_err", condition = is_character(tissue_plot))

          prot_plot <- my_prot_plot(gene_id())
          if (is_character(prot_plot)) {
            output$prot_err <- renderText(prot_plot)
          } else {
            output$prot_plot <- renderPlot(prot_plot)
          }
          toggleElement(
            id = "prot_plot",
            condition = !is_character(prot_plot)
          )
          toggleElement(id = "prot_err", condition = is_character(prot_plot))
        }
      })
    }
  )
}