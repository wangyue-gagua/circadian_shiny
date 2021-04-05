

tableUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    h4("Go information"),
    DT::dataTableOutput(ns('go')),
    h4("homologs in ", em("arabidopsis thaliana")),
    DT::dataTableOutput(ns('homo_ara'))
  )
}

tableServer <- function(id, plot_event, gene_id) {
  moduleServer(id,
               
               function(input, output, session) {
                 # default
                 de_go <- my_go_table("Ghir_D11G029140")
                 
                 output$go <- DT::renderDataTable(de_go)
                 output$homo_ara <-
                   DT::renderDataTable(expr = datatable((
                     my_ara_homo_tbl("Ghir_D11G029140") %>%
                       mutate(Match = map(Match,  ~ as.character(
                         a(
                           href = str_c(
                             "https://www.arabidopsis.org/servlets/TairObject?type=locus&name=",
                             str_match(.x, "(^.*)\\.")[, 2]
                           ),
                           target = "_blank",
                           .x
                         )
                       )))
                   ), escape = FALSE))
                 
                 observeEvent(plot_event(), {
                   if (!str_detect(gene_id(),
                                   "^Ghir_\\w\\d{2}G\\d{6}$|^Ghir_\\w\\d{2}G\\d{6}\\.\\d$")) {
                     showNotification("Invalid input id!", type = "error")
                   } else{
                     go_table <- my_go_table(gene_id())
                     output$go <- DT::renderDataTable(go_table)
                     
                     df <- my_ara_homo_tbl(gene_id())
                     if (!is_character(df))
                       df <- df %>%
                       mutate(Match = map(Match,  ~ as.character(
                         a(
                           href = str_c(
                             "https://www.arabidopsis.org/servlets/TairObject?type=locus&name=",
                             str_match(.x, "(^.*)\\.")[, 2]
                           ),
                           target = "_blank",
                           .x
                         )
                       )))
                     output$homo_ara <-
                       DT::renderDataTable(expr = datatable(df, escape = FALSE))
                   }
                 })
               })
}
