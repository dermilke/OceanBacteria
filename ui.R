### User Interface for OceanBacteria App ###

ui = fluidPage(
  title = "Microbial Abundance Patterns of the Oceans",
  tabsetPanel(
    tabPanel("Browse Dataset",
      sidebarLayout(
        sidebarPanel(
          selectizeInput(inputId = "Phylum",
                         label = "Select Phylum",
                         choices = NULL,
                         options = list(placeholder = 'select')
          ),
          selectizeInput(inputId = "Class",
                         label = "Select Class",
                         choices = NULL,
                         options = list(placeholder = 'select',
                                        onInitialize = I('function() { this.setValue(""); }'))
          ),
          selectizeInput(inputId = "Order",
                         label = "Select Order",
                         choices = NULL,
                         options = list(placeholder = 'select',
                                        onInitialize = I('function() { this.setValue(""); }'))
          ),
          selectizeInput(inputId = "Family",
                         label = "Select Family",
                         choices = NULL,
                         options = list(placeholder = 'select',
                                        onInitialize = I('function() { this.setValue(""); }'))
          ),
          selectizeInput(inputId = "Genus",    
                         label = "Select Genus",
                         choices = NULL,
                         options = list(placeholder = 'select',
                                        onInitialize = I('function() { this.setValue(""); }'))
          ),
          selectizeInput(inputId = "Species",
                         label = "Select Species",
                         choices = NULL,
                         options = list(placeholder = 'select',
                                        onInitialize = I('function() { this.setValue(""); }'))
          ),
          hr(),
          checkboxInput("AbundanceFilter", "Abundance Filter", FALSE)
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Atlantic", plotOutput("AtlanticBrowse")),
            tabPanel("Pacific", plotOutput("PacificBrowse"))
          ),
          hr(),
          DT::dataTableOutput("ASV")
        )
      )
    ),
    tabPanel("BLAST Search",
      sidebarLayout(
        sidebarPanel(
          textInput("seqInputBlast", "Blast Sequence"),
          actionButton("blastGo", "Blast")
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Atlantic", plotOutput("AtlanticBlast")),
            tabPanel("Pacific", plotOutput("PacificBlast"))
          ),
          hr(),
          DT::dataTableOutput("BLAST")
        )
      )
    ),
    tabPanel("Correlation Analysis",
      sidebarLayout(
        sidebarPanel(
          tags$p("This section shows biplots for defined taxa. First, select an OTU from one of the tables, then choose >Abundance< in one of the drop-down-menus."),
          tags$p("Beware priority: If taxa in both tables are selected, only the >Browse Dataset< selection will be displayed under >Abundance<."),
          tags$p(""),
          tags$p("Plotting may take some seconds."),
          hr(),
          panel(heading = "Plot 1",
            varSelectInput("Xparam", "Environment 1 - X value", Meta_Combined_working,
                           selected = "Pot_Temperature"),
            varSelectInput("Yparam", "Environment 2 - Y value", Meta_Combined_working,
                           selected = "Abundance"),
            varSelectInput("Colparam", "Element - Color", Meta_Combined_working,
                           selected = "Province"),
            varSelectInput("Sizeparam", "Element - Size", Meta_Combined_working,
                           selected = "Depth")
          )
        ),
        mainPanel(
          plotOutput("PlotCorAnalysis")
        )
      )
    )
  )
)