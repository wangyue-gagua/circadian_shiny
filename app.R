

library(shiny)
library(tidyverse)
library(readxl)
library(shinythemes)
# library(gridExtra)


# import sample list , raw gene count and TMM normalized expression matrix
# genes.counts <- read.delim("./merged_counts/genes.counts.matrix", row.names=1, check.names = FALSE )
# genes.TMM.EXPR <- read.delim("./merged_counts/genes.TMM.EXPR.matrix", row.names=1, check.names = FALSE)
sample_info <-
  read.table("./merged_counts/sample.info",
             row.names = 'sample',
             header = TRUE)
sample_info_exp <- read_csv("./merged_counts/sample_info_exp.csv")

iso_exp_tpm <- read_delim(
  "data/isoforms(HUA-ccNET).fpkm_table",
  "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)
AD_pro_19920 <- read_excel("data/AD-pro-19920.xlsx")
blastp_AD1_HAU_v1_0_vs_arabidopsis_1 <-
  read_delim(
    "./bla_go/blastp_AD1_HAU_v1.0_vs_arabidopsis.1.txt",
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
genes2Go <- read_excel("./bla_go/genes2Go.xlsx",
                       col_names = T,
                       skip = 1)

# careate an integrated exp and info table
# tran_exp <- as_tibble(t(genes.TMM.EXPR), rownames = 'sample')
### arrange by time
sample_info <-
  sample_info %>% rownames_to_column(var = 'sample') %>% group_by(strain) %>%
  arrange(time, .by_group = TRUE)

trans_time <- function(str, pattr) {
  mat_obj <- str_match(str, pattern = pattr)
  times <- if_else(mat_obj[[3]] == "minus", '-', '')
  day <- if_else(mat_obj[[5]] == "12", 'N', 'L')
  out_str <- str_c(times, mat_obj[[4]], "DPA", '-', day)
  return(out_str)
}
sample_info <-
  sample_info %>% ungroup() %>% mutate(labs = map_chr(.$sample, trans_time, "(\\w+)-(\\w*)(\\d)-DPA(\\d+)h-(\\d)"))

# genes.TMM.EXPR <- genes.TMM.EXPR %>% select(sample_info$sample)

# sample_info_exp <- right_join(sample_info, tran_exp, by='sample' )
# sample_info_exp %>% write_csv("./merged_counts/sample_info_exp.csv")

## reduce replication ,set mean of replication as exp value, a new function
rep_reduce <- function(st) {
  # geneid as input , df with mean value and standard deviation as output
  mean_li <- c()
  for (i in 1:(length(sample_info_exp[[st]]) / 2)) {
    mean_li[i] <-
      mean(c(sample_info_exp[[st]][[2 * i - 1]], sample_info_exp[[st]][[2 * i]]))
  }
  
  std_li <- c()
  for (i in 1:(length(sample_info_exp[[st]]) / 2)) {
    std_li[i] <-
      sd(c(sample_info_exp[[st]][[2 * i - 1]], sample_info_exp[[st]][[2 * i]]))
  }
  # m <- str_c('mean_',st)
  # s <- str_c('std_',st)
  temp <- sample_info %>% filter(replicate == 1) %>% ungroup() %>%
    # mutate(!!m:= mean_li, !!s:= std_li)
    mutate(mean_li, std_li)
  return(temp)
}
## the plot function
my_cir_plot <- function(geneid, alia_name = NULL) {
  # if(alia_name==''){
  #   alia_name <- geneid
  # }
  rects <-
    data.frame(xstart = seq(-36, 36, 24), xend = seq(-24, 48, 24))
  sub_df <- rep_reduce(geneid)
  ggplot(sub_df, aes(x = time, y = mean_li, color = strain)) +
    geom_point(size = 3) +
    # scale_x_continuous(breaks = seq(-72,72,24))+
    facet_grid(rows = vars(strain)) +
    geom_errorbar(
      aes(ymax = std_li + mean_li,
          ymin = mean_li - std_li),
      width = 2,
      color = 'black'
    ) +
    geom_line() +
    ylab("relative expression/(TMM)") +
    scale_x_continuous(breaks = seq(-48, 48, 12), labels = sub_df$labs[1:(length(sub_df$labs) /
                                                                            2)]) +
    geom_rect(
      data = rects,
      aes(
        xmin = xstart,
        xmax = xend,
        ymin = 0,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    labs(title = geneid, subtitle = alia_name)
  
}

my_tissue_explot <- function(st) {
  tryCatch(
    error = function(cnd) {
      str_c("No data available! ", st)
    },
    {
      temp <-
        iso_exp_tpm %>% filter(str_detect(tracking_id, st)) %>% column_to_rownames(var = "tracking_id") %>%
        mutate("leaf_0_" = leaf_0) %>% select(!leaf_0)
      if(nrow(temp)>0){temp <- temp[1,]}
      df <-
        tibble(
          "samples" = factor(names(temp), levels = names(temp)),
          "fpkm" = as.numeric(temp),
          "tissue" = map_chr(names(temp), ~ str_match(.x, "(\\w+?)_.*")[[2]])
        )
      df %>% ggplot(aes(x = samples, y = log(fpkm + 1))) +
        geom_point(aes(color = tissue), size = 3) +
        geom_path(aes(
          group = tissue,
          y = log(fpkm + 1),
          color = tissue
        )) +
        theme(axis.text.x = element_text(angle = 90)) +
        ylab("relative expression/log(fpkm+1)") +
        labs(title = st)
    }
  )
}

my_prot_plot <- function(st) {
  tryCatch(
    error = function(cnd) {
      str_c("No data available! ", st)
    },
    
    {
      df <- AD_pro_19920 %>% filter(str_detect(tracking_id, st))
      if(nrow(df)>0){df <- df[1,]}
      df_pl <-
        tibble(
          samples = factor(
            df %>% select(contains("WT")) %>% names(),
            levels = (df %>% select(contains("WT")) %>% names())
          ),
          exps = (df %>% select(contains("WT")) %>% as_vector())
        )
      ggplot(data = df_pl, aes(
        x = samples,
        y = exps,
        color = str_detect(samples, "tmt")
      )) +
        geom_point(size = 3) +
        labs(
          title = st,
          subtitle = str_c("phosphorylation site: ", df$sites),
          col = "treatment"
        ) +
        theme(legend.text = element_blank())
    }
  )
}

my_go_table <- function(str) {
  str_gtl <- genes2Go %>% filter(str_detect(Query, pattern = str))
  if (nrow(str_gtl) == 0)
    map_dfr(str_gtl, function(x)
      x <- "No GO result")
  else
    return(str_gtl)
}

my_ara_homo_tbl <- function(str) {
  str_homo <-
    blastp_AD1_HAU_v1_0_vs_arabidopsis_1 %>% filter(str_detect(Query, pattern = str))
  if (nrow(str_homo) == 0)
    return("No homolog in tair10")
  else
    return(str_homo)
}

# data preprocess____________________________________________________________


# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("simplex"),
  # Application title
  titlePanel("Circadian rhythm plot generator"),
  

  sidebarLayout(
    sidebarPanel(
      helpText(
        "Input a gene ID, e.g",
        em("Ghir_D11G029140, Ghir_D11G029140.1"),
        br(),
        'click plot for expression, ',
        'GO information and homologs founded in TAIR10 '
      ),
      textInput(
        "geneid",
        label = h3("Input a gene ID"),
        placeholder = "Ghir_D11G029140"
      ),
      actionButton('plot', 'Plot'),
      downloadButton("downloadfig", "Download png"),
      width = 3
    ),
    
    mainPanel(tabsetPanel(
      tabPanel(
        "plots",
        h4("circadian rhythm, gene expression level"),
        plotOutput("circa_plot"),
        h4("expression level of isoform in different tissues"),
        plotOutput("tissue_plot"),
        textOutput("tissue_err", container = h1),
        h4("different expression level of protein between FL and WT"),
        plotOutput("prot_plot"),
        textOutput("prot_err", container = h1)
      ),
      tabPanel(
        "tables",
        h4("Go information"),
        DT::dataTableOutput('go'),
        h4("homologs in ", em("arabidopsis thaliana")),
        DT::dataTableOutput('homo_ara')
      )
    ))
    
    # Show a plot of the generated distribution
  )
)

#Define server logic required to draw a histogram
server <- function(input, output) {
  require(DT)
  
  # default
  output$circa_plot <- renderPlot(my_cir_plot("Ghir_D11G029140"))
  output$tissue_plot <- renderPlot(my_tissue_explot("Ghir_D11G029140"))
  output$prot_plot <- renderPlot(my_prot_plot("Ghir_D11G029140"))
  output$go <- DT::renderDataTable(my_go_table("Ghir_D11G029140"))
  output$homo_ara <-
    DT::renderDataTable(expr = datatable((my_ara_homo_tbl("Ghir_D11G029140")%>%
                                            mutate(Match = map(Match,  ~ as.character(a(
                                              href = str_c(
                                                "https://www.arabidopsis.org/servlets/TairObject?type=locus&name=",
                                                str_match(.x, "(^.*)\\.")[, 2]
                                              ),
                                              target = "_blank",
                                              .x
                                            ))))), escape = FALSE))
  
  observeEvent(input$plot,
               {if(!str_detect(input$geneid, "^Ghir_\\w\\d{2}G\\d{6}$|^Ghir_\\w\\d{2}G\\d{6}\\.\\d$")){
                 showNotification("Invalid input id!", type = "error")
               }else{
                 cir_plot <- my_cir_plot(str_extract(input$geneid, "^Ghir_\\w\\d{2}G\\d{6}"))
                 output$circa_plot <- renderPlot(cir_plot)
                 
                 tissue_plot <- my_tissue_explot(input$geneid)
                 if (is_character(tissue_plot)) {
                   output$tissue_err <- renderText(tissue_plot)
                   output$tissue_plot <- renderPlot(NULL)
                 } else{
                   output$tissue_plot <- renderPlot(tissue_plot)
                   output$tissue_err <- renderText(NULL)
                 }
                 
                 prot_plot <- my_prot_plot(input$geneid)
                 if (is_character(prot_plot)) {
                   output$prot_err <- renderText(prot_plot)
                   output$prot_plot <- renderPlot(NULL)
                 } else{
                   output$prot_plot <- renderPlot(prot_plot)
                   output$prot_err <- renderText(NULL)
                 }
                 
                 go_table <- my_go_table(input$geneid)
                 output$go <- DT::renderDataTable(go_table)
                 
                 df <- my_ara_homo_tbl(input$geneid)
                 if (!is.null(df$Match))
                   df <- df %>%
                   mutate(Match = map(Match,  ~ as.character(a(
                     href = str_c(
                       "https://www.arabidopsis.org/servlets/TairObject?type=locus&name=",
                       str_match(.x, "(^.*)\\.")[, 2]
                     ),
                     target = "_blank",
                     .x
                   ))))
                 output$homo_ara <-
                   DT::renderDataTable(expr = datatable(df, escape = FALSE))
               }})
  
  # Downloadable csv of selected dataset ----
  output$downloadfig <- downloadHandler(
    filename = function() {
      paste(input$geneid, ".png", sep = "")
    },
    content = function(file) {
      ggsave(filename = file,
             plot = my_cir_plot(input$geneid),
             dpi = 500)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
