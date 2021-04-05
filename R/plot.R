

# data preprocess____________________________________________________________



#_____________________________________________________________________________

plotUI <- function(id) {
  ns <- NS(id)
  tagList(
      h4("circadian rhythm, gene expression level"),
      plotOutput(ns("circa_plot")),
      textOutput(ns('circa_err'), container = h1),
      h4("expression level of isoform in different tissues"),
      plotOutput(ns("tissue_plot")),
      textOutput(ns("tissue_err"), container = h1),
      h4("different expression level of protein between FL and WT"),
      plotOutput(ns("prot_plot")),
      textOutput(ns("prot_err"), container = h1)
  )
  
}

plotServer <- function(id, plot_event, gene_id) {
  moduleServer(id,
               
               function(input, output, session) {
                 ###
                 # default plot page
                 de_cir <- my_cir_plot("Ghir_D11G029140")
                 
                 de_tissue <- my_tissue_plot("Ghir_D11G029140")
                 de_pro <- my_prot_plot("Ghir_D11G029140")
                 
                 
                 output$circa_plot <- renderPlot(de_cir)
                 output$tissue_plot <-
                   renderPlot(de_tissue)
                 output$prot_plot <- renderPlot(de_pro)

                 observeEvent(plot_event(), {
                   if (!str_detect(gene_id(),
                                   "^Ghir_\\w\\d{2}G\\d{6}$|^Ghir_\\w\\d{2}G\\d{6}\\.\\d$")) {
                     showNotification("Invalid input id!", type = "error")
                   } else{
                     circa_plot <-
                       my_cir_plot(str_extract(gene_id(), "^Ghir_\\w\\d{2}G\\d{6}"))
                     
                     
                     if (is_character(circa_plot)) {
                       output$circa_err <- renderText(circa_plot)
                       output$circa_plot <- renderPlot(NULL)
                       hide("circa_plot")
                     } else{
                       output$circa_plot <- renderPlot(circa_plot)
                       output$circa_err <- renderText(NULL)
                       show("circa_plot")
                     }
                     
                     tissue_plot <- my_tissue_plot(gene_id())
                     if (is_character(tissue_plot)) {
                       output$tissue_err <- renderText(tissue_plot)
                       output$tissue_plot <- renderPlot(NULL)
                       hide("tissue_plot")
                     } else{
                       output$tissue_plot <- renderPlot(tissue_plot)
                       output$tissue_err <- renderText(NULL)
                       show("tissue_plot")
                     }
                     
                     
                     prot_plot <- my_prot_plot(gene_id())
                     if (is_character(prot_plot)) {
                       output$prot_err <- renderText(prot_plot)
                       output$prot_plot <- renderPlot(NULL)
                       hide("prot_plot")
                     } else{
                       output$prot_plot <- renderPlot(prot_plot)
                       output$prot_err <- renderText(NULL)
                       show("prot_plot")
                     }
                     
                     
                   }
                 })
                 
  
               })
}
