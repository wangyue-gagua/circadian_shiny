



library(shiny)
library(shinyjs)
library(tidyverse)
library(readxl)
library(shinythemes)
library(svglite)
library(pheatmap)
# library(gridExtra)


# import sample list , raw gene count and TMM normalized expression matrix
# genes.counts <- read.delim("./merged_counts/genes.counts.matrix", row.names=1, check.names = FALSE )
genes.TMM.EXPR <-
  read.delim(
    "./merged_counts/genes.TMM.EXPR.matrix",
    row.names = 1,
    check.names = FALSE
  )

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
#### trans time to labels
trans_time <- function(str, pattr) {
  mat_obj <- str_match(str, pattern = pattr)
  times <- if_else(mat_obj[[3]] == "minus", '-', '')
  day <- if_else(mat_obj[[5]] == "12", 'N', 'L')
  out_str <- str_c(times, mat_obj[[4]], "DPA", '-', day)
  return(out_str)
}
sample_info <-
  sample_info %>% ungroup() %>% mutate(labs = map_chr(.$sample, trans_time, "(\\w+)-(\\w*)(\\d)-DPA(\\d+)h-(\\d)"))

genes.TMM.EXPR <- genes.TMM.EXPR %>% select(sample_info$sample)

# sample_info_exp <- right_join(sample_info, tran_exp, by='sample' )
# sample_info_exp %>% write_csv("./merged_counts/sample_info_exp.csv")

## reduce replication ,set mean of replication as exp value, a new function
rep_reduce <- function(st) {
  tempvar <- sym(st)
  sample_info_exp %>% select(c(1:6), !!tempvar,) %>% group_by(strain, period, time, labs) %>%
    summarise(
      'mean_li' = mean(!!tempvar),
      'std_li' = sd(!!tempvar),
      .groups = "keep"
    )
}
## the plot function
my_cir_plot <- function(geneid, alia_name = NULL) {
  tryCatch(
    error = function(cnd) {
      str_c("No RNA-seq data available! ", geneid)
    },
    {
      rects <-
        data.frame(xstart = seq(-36, 36, 24),
                   xend = seq(-24, 48, 24))
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
    
  )
}


my_tissue_plot <- function(st) {
  tryCatch(
    error = function(cnd) {
      str_c("No RNA-seq data available! ", st)
    },
    {
      temp <-
        iso_exp_tpm %>% filter(str_detect(tracking_id, st)) %>% column_to_rownames(var = "tracking_id") %>%
        mutate("leaf_0_" = leaf_0) %>% select(!leaf_0)
      if (nrow(temp) > 0) {
        temp <- temp[1,]
      } else{
        stop("No data available!")
      }
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
      str_c("No MS data available! ", st)
    },
    
    {
      df <- AD_pro_19920 %>% filter(str_detect(tracking_id, st))
      if (nrow(df) > 0) {
        df <- df[1,]
      }
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
    map_dfr(str_homo, function(x)
      x <- "No homolog in tair10")
  else
    return(str_homo)
}


# data preprocess____________________________________________________________

###
# default plot page
de_cir <- my_cir_plot("Ghir_D11G029140")

de_tissue <- my_tissue_plot("Ghir_D11G029140")
de_pro <- my_prot_plot("Ghir_D11G029140")
de_go <- my_go_table("Ghir_D11G029140")
# Define UI for application that draws a histogram
ui <- fluidPage(tabsetPanel(
  tabPanel(
    "generator",
    useShinyjs(),
    # Include shinyjs
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
          'GO information and homologs found in TAIR10 '
        ),
        textInput(
          "geneid",
          label = h3("Input a gene ID"),
          placeholder = "Ghir_D11G029140"
        ),
        #button
        actionButton('plot', 'Plot'),
        
        downloadButton("downloadSourceData", "Download data"),
        downloadButton("downloadfig", "Download fig"),
        selectInput(
          "fileType",
          label = "picture format",
          choices = list(
            "png" = ".png",
            "jpeg" = ".jpeg",
            "svg" = ".svg",
            "pdf" = ".pdf"
          )
        ),
        fileInput('file_load', 'gene id in tab-delimited format', accept = 'text/plain'),
        ### batch tab is active
        uiOutput('pltHtMap'),
        width = 3
      ),
      
      mainPanel(
        tabsetPanel(
          tabPanel(
            title = "plots",
            h4("circadian rhythm, gene expression level"),
            plotOutput("circa_plot"),
            textOutput('circa_err', container = h1),
            h4("expression level of isoform in different tissues"),
            plotOutput("tissue_plot"),
            textOutput("tissue_err", container = h1),
            h4("different expression level of protein between FL and WT"),
            plotOutput("prot_plot"),
            textOutput("prot_err", container = h1)
          ),
          tabPanel(
            title = "tables",
            h4("Go information"),
            DT::dataTableOutput('go'),
            h4("homologs in ", em("arabidopsis thaliana")),
            DT::dataTableOutput('homo_ara')
          ),
          tabPanel(
            title = 'batch',
            h4('heatmap for multiple genes'),
            plotOutput('heatmap'),
          ),
          id = 'main_tabsets'
        )
      )
      # Show a plot of the generated distribution
    )
  ),
  tabPanel(
    "about",
    h1("any problems please contact to wangyue"),
    br(),
    h2("goto github for change log"),
    h3(
      "github: ",
      a(href = "https://github.com/wangyue-gagua/circadian_shiny", "https://github.com/wangyue-gagua/circadian_shiny")
    ),
    br(),
    h3(str_c("last update time: ", date()))
  )
))



#Define server logic required to draw a histogram
server <- function(input, output) {
  require(DT)
  
  # default
  output$circa_plot <- renderPlot(de_cir)
  output$tissue_plot <-
    renderPlot(de_tissue)
  output$prot_plot <- renderPlot(de_pro)
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
  observeEvent(input$plot, {
    if (!str_detect(input$geneid,
                    "^Ghir_\\w\\d{2}G\\d{6}$|^Ghir_\\w\\d{2}G\\d{6}\\.\\d$")) {
      showNotification("Invalid input id!", type = "error")
    } else{
      circa_plot <-
        my_cir_plot(str_extract(input$geneid, "^Ghir_\\w\\d{2}G\\d{6}"))
      
      if (is_character(circa_plot)) {
        output$circa_err <- renderText(circa_plot)
        output$circa_plot <- renderPlot(NULL)
        hide("circa_plot")
      } else{
        output$circa_plot <- renderPlot(circa_plot)
        output$circa_err <- renderText(NULL)
        show("circa_plot")
      }
      
      tissue_plot <- my_tissue_plot(input$geneid)
      if (is_character(tissue_plot)) {
        output$tissue_err <- renderText(tissue_plot)
        output$tissue_plot <- renderPlot(NULL)
        hide("tissue_plot")
      } else{
        output$tissue_plot <- renderPlot(tissue_plot)
        output$tissue_err <- renderText(NULL)
        show("tissue_plot")
      }
      
      
      prot_plot <- my_prot_plot(input$geneid)
      if (is_character(prot_plot)) {
        output$prot_err <- renderText(prot_plot)
        output$prot_plot <- renderPlot(NULL)
        hide("prot_plot")
      } else{
        output$prot_plot <- renderPlot(prot_plot)
        output$prot_err <- renderText(NULL)
        show("prot_plot")
      }
      
      go_table <- my_go_table(input$geneid)
      output$go <- DT::renderDataTable(go_table)
      
      df <- my_ara_homo_tbl(input$geneid)
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
  
  observeEvent(input$file_load,{
    output$pltHtMap <- renderUI({
      actionButton('heat_plot', 'HeatMap')})
  })
 

  
  
  output$heatmap <- renderPlot({
    id_list <-
      read_lines(input$file_load$datapath) %>% str_split('\t') %>% as_vector()
    if_show_rownames = if_else(length(id_list) < 20, true = TRUE, false = FALSE)
    df <- genes.TMM.EXPR[id_list, ]
    pheatmap(
      df %>% na.omit(),
      cluster_cols = FALSE,
      show_rownames = if_show_rownames,
      annotation_col = sample_info %>% column_to_rownames(var =
                                                            'sample'),
      scale = 'none'
    )
  }) %>%
    bindCache(input$file_load$datapath) %>%
    bindEvent(input$heat_plot)
  
  
  
  
  
  
  # Downloadable csv of selected dataset
  output$downloadfig <- downloadHandler(
    filename = function() {
      paste(input$geneid, input$fileType, sep = "")
    },
    content = function(file) {
      ggsave(filename = file,
             plot = my_cir_plot(input$geneid))
    }
  )
  
  output$downloadSourceData <- downloadHandler(
    filename = function() {
      paste(input$geneid, ".csv", sep = "")
    },
    content = function(filename) {
      sample_info_exp %>% select(sample, strain, period, time, replicate, labs, input$geneid) %>%
        write_csv(file = filename, col_names = TRUE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
