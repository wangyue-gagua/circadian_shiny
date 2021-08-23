source("./helper/data_load.R")
source("./helper/method.R")


ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "singleCell.css")
  ),
  tabsetPanel(
  tabPanel(
    "generator",
    useShinyjs(),
    # Include shinyjs
    theme = shinytheme("simplex"),
    # Application title
    titlePanel("Circadian rhythm plot generator"),
    sidebarLayout(
      sidebarPanel(
        geneInputUI("input_id"),
        downloadUI("download"),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            title = "plots",
            plotUI("plots")
          ),
          tabPanel(
            title = "tables",
            tableUI("tables")
          ),
          tabPanel(
            title = "single_cell",
            singleCellUI("single_cell")
          ),
          tabPanel(
            title = "batch",
            batchUI("batch")
          )
        )
      )
      # Show a plot of the generated distribution
    )
  ),
  tabPanel(
    "tools",
    sidebarLayout(
      sidebarPanel(
        IdConvert_UI("ConvertIn"),
        width = 3
      ),
      mainPanel(
        ConvertOut_UI("ConvertOut")
      )
    )
  ),
  tabPanel(
    "about",
    h1("any problems please contact to wangyue"),
    br(),
    h2("goto github for change log"),
    h3(
      "github: ",
      a(href = "https://github.com/wangyue-gagua/circadian_shiny",
       "https://github.com/wangyue-gagua/circadian_shiny")
    ),
    br(),
    h3(str_c("last update time: ", date()))
  )
))



# Define server logic required to draw a histogram
server <- function(input, output) {
  require(DT)

  input_id <- geneInputServer("input_id")
  plotServer("plots", input_id$tri, input_id$id)
  tableServer("tables", input_id$tri, input_id$id)
  singleCellServer("single_cell", input_id$tri, input_id$id)

  batchServer("batch")

  downloadServer("download", input_id$id)
  ## tools
  convertOpt <- IdConvert_Server("ConvertIn")
  ConvertOutServer("ConvertOut", convertOpt$genome ,convertOpt$id, convertOpt$run)
}

# Run the application
shinyApp(ui = ui, server = server)