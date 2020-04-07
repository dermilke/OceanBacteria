### User Interface for OceanBacteria App ###

ui = fluidPage(
  title = "Microbial Abundance Patterns of the Oceans",
  tabsetPanel(
    #### Prokaryotes Panel ####
    tabPanel("Prokaryotes",
      tabsetPanel(
        #### Browse Prokaryotes ####
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
                                            onInitialize = I('function() { this.setValue("All"); }'))
              ),
              selectizeInput(inputId = "Order",
                             label = "Select Order",
                             choices = NULL,
                             options = list(placeholder = 'select',
                                            onInitialize = I('function() { this.setValue("All"); }'))
              ),
              selectizeInput(inputId = "Family",
                             label = "Select Family",
                             choices = NULL,
                             options = list(placeholder = 'select',
                                            onInitialize = I('function() { this.setValue("All"); }'))
              ),
              selectizeInput(inputId = "Genus",    
                             label = "Select Genus",
                             choices = NULL,
                             options = list(placeholder = 'select',
                                            onInitialize = I('function() { this.setValue("All"); }'))
              ),
              selectizeInput(inputId = "Species",
                             label = "Select Species",
                             choices = NULL,
                             options = list(placeholder = 'select',
                                            onInitialize = I('function() { this.setValue("All"); }'))
              ),
              hr(),
              checkboxInput("AbundanceFilter", "Abundance Filter", FALSE)
            ),
            mainPanel(
              tabsetPanel(
                tabPanel("Map", 
                         fluidRow(
                           column(4, 
                              selectInput(inputId = "MapDepth", label = "Select Depth", 
                                          choices = c(20, 40, 60, 100, 200, 300, 500))),
                           column(4,
                             selectInput(inputId = "MapSizeFraction", label = "Select Size Fraction", 
                                         choices = c("0.22-3µm (FL)", "3-8µm (SPA)", ">8µm (LPA)"))),
                           column(4,
                             sliderInput(inputId = "ColorRange", label = "Define Color Range",
                                         min = 0, max = 0.5, value = c(0,0.2), step = 0.001))),
                         leaflet::leafletOutput("MapOutput")),
                tabPanel("Atlantic", plotOutput("AtlanticBrowse")),
                tabPanel("Pacific", plotOutput("PacificBrowse"))
              ),
              hr(),
              DT::dataTableOutput("ASV")
            )
          )
        ),
        #### BLAST Prokaryotes ####
        tabPanel("BLAST Search",
          sidebarLayout(
            sidebarPanel(
              textInput("seqInputBlast", "Blast Sequence"),
              actionButton("blastGo", "Blast"),
              hr(),
              p("Enter OTU ID here to get the corresponding sequence:"),
              textInput("hashInput", "OTU Identifier"),
              actionButton("invBlastGo", "Get Sequence"),
              verbatimTextOutput("hashOutput")
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
        #### Correlation Prokaryotes ####
        tabPanel("Correlation Analysis",
          sidebarLayout(
            sidebarPanel(
              tags$p("This section shows biplots for defined taxa. First, select an OTU from one of the tables, then choose >Abundance< in one of the drop-down-menus."),
              tags$p("Beware priority: If taxa in both tables are selected, only the >Browse Dataset< selection will be displayed under >Abundance<."),
              tags$p(""),
              tags$p("Plotting may take some seconds."),
              hr(),
              panel(heading = "Plot 1",
                varSelectInput("Xparam", "Environment 1 - X value", datalist_Combined$Meta_Data,
                               selected = "Pot_Temperature"),
                varSelectInput("Yparam", "Environment 2 - Y value", datalist_Combined$Meta_Data,
                               selected = "Abundance"),
                varSelectInput("Colparam", "Element - Color", datalist_Combined$Meta_Data,
                               selected = "Cruise"),
                varSelectInput("Sizeparam", "Element - Size", datalist_Combined$Meta_Data,
                               selected = "Depth")
              )
            ),
            mainPanel(
              plotOutput("PlotCorAnalysis")
            )
          )
        )
      )
    ),
    #### Chloroplasts Panel ####
    tabPanel("Chloroplasts",
             tabsetPanel(
               #### Browse Chloroplasts ####
               tabPanel("Browse Dataset",
                        sidebarLayout(
                          sidebarPanel(
                            selectizeInput(inputId = "ChloroSupergroup",
                                           label = "Select Supergroup",
                                           choices = NULL,
                                           options = list(placeholder = 'select')
                            ),
                            selectizeInput(inputId = "ChloroPhylum",
                                           label = "Select Phylum",
                                           choices = NULL,
                                           options = list(placeholder = 'select',
                                                          onInitialize = I('function() { this.setValue("All"); }'))
                            ),
                            selectizeInput(inputId = "ChloroClass",
                                           label = "Select Class",
                                           choices = NULL,
                                           options = list(placeholder = 'select',
                                                          onInitialize = I('function() { this.setValue("All"); }'))
                            ),
                            selectizeInput(inputId = "ChloroSubclass",
                                           label = "Select Subclass",
                                           choices = NULL,
                                           options = list(placeholder = 'select',
                                                          onInitialize = I('function() { this.setValue("All"); }'))
                            ),
                            selectizeInput(inputId = "ChloroOrder",    
                                           label = "Select Order",
                                           choices = NULL,
                                           options = list(placeholder = 'select',
                                                          onInitialize = I('function() { this.setValue("All"); }'))
                            ),
                            selectizeInput(inputId = "ChloroSuborder",
                                           label = "Select Suborder",
                                           choices = NULL,
                                           options = list(placeholder = 'select',
                                                          onInitialize = I('function() { this.setValue("All"); }'))
                            ),
                            selectizeInput(inputId = "ChloroFamily",
                                           label = "Select Family",
                                           choices = NULL,
                                           options = list(placeholder = 'select',
                                                          onInitialize = I('function() { this.setValue("All"); }'))
                            ),
                            selectizeInput(inputId = "ChloroGenus",
                                           label = "Select Genus",
                                           choices = NULL,
                                           options = list(placeholder = 'select',
                                                          onInitialize = I('function() { this.setValue("All"); }'))
                            ),
                            selectizeInput(inputId = "ChloroSpecies",
                                           label = "Select Species",
                                           choices = NULL,
                                           options = list(placeholder = 'select',
                                                          onInitialize = I('function() { this.setValue("All"); }'))
                            ),
                            hr(),
                            checkboxInput("ChloroAbundanceFilter", "Abundance Filter", FALSE)
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Map", 
                                       fluidRow(
                                         column(4, 
                                                selectInput(inputId = "ChloroMapDepth", label = "Select Depth", 
                                                            choices = c(20, 40, 60, 100, 200, 300, 500))),
                                         column(4,
                                                selectInput(inputId = "ChloroMapSizeFraction", label = "Select Size Fraction", 
                                                            choices = c("0.22-3µm (FL)", "3-8µm (SPA)", ">8µm (LPA)"))),
                                         column(4,
                                                sliderInput(inputId = "ChloroColorRange", label = "Define Color Range",
                                                            min = 0, max = 0.5, value = c(0,0.2), step = 0.001))),
                                       leaflet::leafletOutput("ChloroMapOutput")),
                              tabPanel("Atlantic", plotOutput("ChloroAtlanticBrowse")),
                              tabPanel("Pacific", plotOutput("ChloroPacificBrowse"))
                            ),
                            hr(),
                            DT::dataTableOutput("ChloroASV")
                          )
                        )
               ),
               #### BLAST Chloroplasts ####
               tabPanel("BLAST Search",
                        sidebarLayout(
                          sidebarPanel(
                            textInput("ChloroSeqInputBlast", "Blast Sequence"),
                            actionButton("ChloroBlastGo", "Blast"),
                            hr(),
                            p("Enter OTU ID here to get the corresponding sequence:"),
                            textInput("ChloroHashInput", "OTU Identifier"),
                            actionButton("ChloroInvBlastGo", "Get Sequence"),
                            verbatimTextOutput("ChloroHashOutput")
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Atlantic", plotOutput("ChloroAtlanticBlast")),
                              tabPanel("Pacific", plotOutput("ChloroPacificBlast"))
                            ),
                            hr(),
                            DT::dataTableOutput("ChloroBLAST")
                          )
                        )
               ),
               #### Correlation Chloroplasts ####
               tabPanel("Correlation Analysis",
                        sidebarLayout(
                          sidebarPanel(
                            tags$p("This section shows biplots for defined taxa. First, select an OTU from one of the tables, then choose >Abundance< in one of the drop-down-menus."),
                            tags$p("Beware priority: If taxa in both tables are selected, only the >Browse Dataset< selection will be displayed under >Abundance<."),
                            tags$p(""),
                            tags$p("Plotting may take some seconds."),
                            hr(),
                            panel(heading = "Plot 1",
                                  varSelectInput("ChloroXparam", "Environment 1 - X value", datalist_Combined$Meta_Data,
                                                 selected = "Pot_Temperature"),
                                  varSelectInput("ChloroYparam", "Environment 2 - Y value", datalist_Combined$Meta_Data,
                                                 selected = "Abundance"),
                                  varSelectInput("ChloroColparam", "Element - Color", datalist_Combined$Meta_Data,
                                                 selected = "Cruise"),
                                  varSelectInput("ChloroSizeparam", "Element - Size", datalist_Combined$Meta_Data,
                                                 selected = "Depth")
                            )
                          ),
                          mainPanel(
                            plotOutput("ChloroPlotCorAnalysis")
                          )
                        )
               )
             )
    )
  )
)
