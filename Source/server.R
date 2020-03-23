### Serverside output rendering ###

server = function(input, output, session) {
  
  output$PlotCorAnalysis <- renderPlot({
    
    if (input$Colparam == "Empty" | is_empty(input$Colparam) | input$Colparam == "") {
      colorValue <- "red"
    } else {
      colorValue <- input$Colparam
    }
    
    if (input$Sizeparam == "Empty" | is_empty(input$Sizeparam) | input$Sizeparam == "") {
      sizeValue <- 10
    } else {
      sizeValue <- input$Sizeparam
    }
    
    Meta_Combined_working() %>%
      ggplot(., aes(x = !!input$Xparam, y = !!input$Yparam, col = !!input$Colparam, size = !!input$Sizeparam)) +
      geom_point()
  
  })
  
  Meta_Combined_working <- reactive({
    
    if (is.numeric(input$ASV_rows_selected)) {
      
      datalist_Combined <- combineOceanDatalists(
                             datalist_Atlantic_working()  %>%
                               makeAppProportion(., apply(select_if(datalist_Atlantic$Count_Data, is.numeric), 2, sum)), 
                             datalist_Pacific_working()  %>%
                               makeAppProportion(., apply(select_if(datalist_Pacific$Count_Data, is.numeric), 2, sum))
                           )
      Meta_Combined_working <- datalist_Combined$Meta_Data %>%
                                 mutate(Abundance = datalist_Combined$Count_Data %>%
                                          slice(input$ASV_rows_selected) %>%
                                          select_if(is.numeric) %>%
                                          as.numeric()
                                 )
      
    } else if (is.numeric(input$BLAST_rows_selected)) {
      
      datalist_Combined <- combineOceanDatalists(
                             datalist_Atlantic %>%
                                makeAppProportion(., apply(select_if(datalist_Atlantic$Count_Data, is.numeric), 2, sum)), 
                             datalist_Pacific %>%
                                makeAppProportion(., apply(select_if(datalist_Pacific$Count_Data, is.numeric), 2, sum))
                           )
      Meta_Combined_working <- datalist_Combined$Meta_Data %>%
                                  mutate(Abundance = datalist_Combined$Count_Data %>%
                                    filter(
                                       `#OTU_ID` == {datalist_Blast_working() %>% 
                                       slice(input$BLAST_rows_selected) %>%
                                       select(1) %>%
                                       as_vector()
                                    }) %>%
                                    select_if(is.numeric) %>%
                                    as.numeric()
                                  )
      
    } else {
      datalist_Combined$Meta_Data
    }
  })
  
  datalist_Atlantic_working <- reactive({
    
    tmp <- datalist_Atlantic
    
    if (input$AbundanceFilter) {
      tmp <- tmp %>%
        filterMiliciStyle(.)
    }
    
    if (input$Phylum == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        group_by(Phylum) %>%
        summarize_if(is.numeric, sum)
      
    } else if (input$Class == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        group_by(Class) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$Order == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        filter(Class == input$Class) %>%
        group_by(Order) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$Family == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        filter(Class == input$Class) %>%
        filter(Order == input$Order) %>%
        group_by(Family) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$Genus == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        filter(Class == input$Class) %>%
        filter(Order == input$Order) %>%
        filter(Family == input$Family) %>%
        group_by(Genus) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$Species == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        filter(Class == input$Class) %>%
        filter(Order == input$Order) %>%
        filter(Family == input$Family) %>%
        filter(Genus == input$Genus) %>%
        group_by(Species) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>% 
        filter(Class == input$Class) %>% 
        filter(Order == input$Order) %>% 
        filter(Family == input$Family) %>% 
        filter(Genus == input$Genus) %>%
        filter(Species == input$Species)
      
    }
    tmp
  })
  
  datalist_Pacific_working <- reactive({
    
    tmp <- datalist_Pacific
    
    if (input$AbundanceFilter) {
      tmp <- tmp %>%
        filterMiliciStyle(.)
    }
    
    if (input$Phylum == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        group_by(Phylum) %>%
        summarize_if(is.numeric, sum)
      
    } else if (input$Class == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        group_by(Class) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$Order == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        filter(Class == input$Class) %>%
        group_by(Order) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$Family == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        filter(Class == input$Class) %>%
        filter(Order == input$Order) %>%
        group_by(Family) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$Genus == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        filter(Class == input$Class) %>%
        filter(Order == input$Order) %>%
        filter(Family == input$Family) %>%
        group_by(Genus) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$Species == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>%
        filter(Class == input$Class) %>%
        filter(Order == input$Order) %>%
        filter(Family == input$Family) %>%
        filter(Genus == input$Genus) %>%
        group_by(Species) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Phylum == input$Phylum) %>% 
        filter(Class == input$Class) %>% 
        filter(Order == input$Order) %>% 
        filter(Family == input$Family) %>% 
        filter(Genus == input$Genus) %>%
        filter(Species == input$Species)
      
    }
    
    tmp 
  })
  
  datalist_Taxonomy_working <- reactive({ 
    
    updateTaxonomy(datalist_Atlantic_working(), datalist_Pacific_working())
    
  })
  
  datalist_Blast_working <- eventReactive(input$blastGo, {
    
    if (is.character(input$seqInputBlast) & nchar(input$seqInputBlast) > 4 &
        grepl("^[ATCG]+$", input$seqInputBlast)) {
      
      write_file(paste(">Query\n", input$seqInputBlast, sep = ""), queryFile)
      
      blastTable <- blastWrapper(blastCommand, blastDB, queryFile) %>%
        filter(pident > 0.95)
      
      output <- right_join(datalist_Taxonomy, blastTable, by = "#OTU_ID") %>%
        drop_na()
      
      output
    
    } else {
      
      showNotification("No sequence input for blast or sequence is too short (min length = 4).")
      tibble()
    
    }
  })
  
  updateSelectizeInput(
    session = session,
    inputId = "Phylum",
    choices = c("All", datalist_Taxonomy %>%
                  select(Phylum) %>%
                  as_vector() %>%
                  unique() %>%
                  sort()
    ),
    selected = "",
    options = list(placeholder = 'select'),
    server = TRUE
  )
  
  observeEvent(input$Phylum, {
    choice_Class <- datalist_Taxonomy %>%
      filter(if (input$Phylum != "All") {.$Phylum == input$Phylum} else {TRUE}) %>%
      select(Class) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "Class",
      choices = c("All", choice_Class))
    
  })
  
  observeEvent(input$Class, {
    choice_Order <- datalist_Taxonomy %>%
      filter(if (input$Phylum != "All") {.$Phylum == input$Phylum} else {TRUE}) %>%
      filter(if (input$Class != "All") {.$Class == input$Class} else {TRUE}) %>%
      select(Order) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "Order",
      choices = c("All", choice_Order)
    )
  })
  
  observeEvent(input$Order, {
    choice_Family <- datalist_Taxonomy %>%
      filter(if (input$Phylum != "All") {.$Phylum == input$Phylum} else {TRUE}) %>%
      filter(if (input$Class != "All") {.$Class == input$Class} else {TRUE}) %>%
      filter(if (input$Order != "All") {.$Order == input$Order} else {TRUE}) %>%
      select(Family) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "Family",
      choices = c("All", choice_Family)
    )
  })
  
  observeEvent(input$Family, {
    choice_Genus <- datalist_Taxonomy %>%
      filter(if (input$Phylum != "All") {.$Phylum == input$Phylum} else {TRUE}) %>%
      filter(if (input$Class != "All") {.$Class == input$Class} else {TRUE}) %>%
      filter(if (input$Order != "All") {.$Order == input$Order} else {TRUE}) %>%
      filter(if (input$Family != "All") {.$Family == input$Family} else {TRUE}) %>%
      select(Genus) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "Genus",
      choices = c("All", choice_Genus)
    )
  })
  
  observeEvent(input$Genus,{
    choice_Species <- datalist_Taxonomy %>%
      filter(if (input$Phylum != "All") {.$Phylum == input$Phylum} else {TRUE}) %>%
      filter(if (input$Class != "All") {.$Class == input$Class} else {TRUE}) %>%
      filter(if (input$Order != "All") {.$Order == input$Order} else {TRUE}) %>%
      filter(if (input$Family != "All") {.$Family == input$Family} else {TRUE}) %>%
      select(Species) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "Species",
      choices = c("All", choice_Species)
    )
  })
  
  output$ASV <- DT::renderDataTable(
    DT::datatable(datalist_Taxonomy_working() %>%
                    rename("Counts Atlantic" = "Counts_1",
                           "Counts Pacific" = "Counts_2",
                           "Richness Atlantic" = "Richness_1",
                           "Richness Pacific" = "Richness_2"),
                  selection = "single")
  )
  
  output$BLAST <- DT::renderDataTable(
    DT::datatable(datalist_Blast_working() %>%
                    rename("Counts Atlantic" = "Counts_1",
                           "Counts Pacific" = "Counts_2",
                           "Richness Atlantic" = "Richness_1",
                           "Richness Pacific" = "Richness_2"),
                  selection = "single")
  )
  
  output$PacificBrowse <- renderPlot({
    
    if (is.numeric(input$ASV_rows_selected)) {
      
      if (select(slice(datalist_Taxonomy_working(), input$ASV_rows_selected), Counts_2) > 0) {
      
        inputGrp <- datalist_Taxonomy_working() %>%
          slice(input$ASV_rows_selected) %>%
          select(1) %>%
          as_vector()
        
        ASVname <- paste(ifelse(input$Phylum == "All", inputGrp, input$Phylum),
                         ifelse(input$Class == "All", inputGrp, input$Class),
                         ifelse(input$Order == "All", inputGrp, input$Order),
                         ifelse(input$Family == "All", inputGrp, input$Family),
                         ifelse(input$Genus == "All", inputGrp, input$Genus),
                         ifelse(input$Species == "All", inputGrp, input$Species), sep = " - ")
        
        plotSingle(datalist_Pacific_working() %>%
                     makeAppProportion(., Ocean = "Pacific"), inputGrp, ASVname)
        
      } else {
        
        ggplot(tibble(text = "ASV not found in Pacific.")) +
          geom_text(aes(x = 1, y = 1, label = text)) +
          theme_void()
        
      }
    }
  })
  
  output$AtlanticBrowse <- renderPlot({
    
    if (is.numeric(input$ASV_rows_selected)) {
      
      if (select(slice(datalist_Taxonomy_working(), input$ASV_rows_selected), Counts_1) > 0) {
        
        inputGrp <- datalist_Taxonomy_working() %>%
          slice(input$ASV_rows_selected) %>%
          select(1) %>%
          as_vector()
        
        ASVname <- paste(ifelse(input$Phylum == "All", inputGrp, input$Phylum),
                         ifelse(input$Class == "All", inputGrp, input$Class),
                         ifelse(input$Order == "All", inputGrp, input$Order),
                         ifelse(input$Family == "All", inputGrp, input$Family),
                         ifelse(input$Genus == "All", inputGrp, input$Genus),
                         ifelse(input$Species == "All", inputGrp, input$Species), sep = " - ")
        
        plotSingle(datalist_Atlantic_working() %>%
                     makeAppProportion(., Ocean = "Atlantic"), inputGrp, ASVname)
        
      } else {
        
        ggplot(tibble(text = "ASV not found in Atlantic")) +
          geom_text(aes(x = 1, y = 1, label = text)) +
          theme_void()
        
      }
    }
  })
  
  output$PacificBlast <- renderPlot({
    
    if (is.numeric(input$BLAST_rows_selected)) {
      
      if (select(slice(datalist_Blast_working(), input$BLAST_rows_selected), Counts_2) > 0) {
        
        inputGrp <- datalist_Blast_working() %>%
          slice(input$BLAST_rows_selected)
        
        ASVname <- paste(inputGrp$Phylum, inputGrp$Class, inputGrp$Order,
                         inputGrp$Family, inputGrp$Genus, inputGrp$Species, sep = " - ")
        
        plotSingle(datalist_Pacific %>%
                     makeAppProportion(., Ocean = "Pacific"), as_vector(select(inputGrp,1)), ASVname)
        
      } else {
        
        ggplot(tibble(text = "ASV not found in Pacific.")) +
          geom_text(aes(x = 1, y = 1, label = text)) +
          theme_void()
        
      }
    }
  })
  
  output$AtlanticBlast <- renderPlot({
    
    if (is.numeric(input$BLAST_rows_selected)) {
      
      if (select(slice(datalist_Blast_working(), input$BLAST_rows_selected), Counts_1) > 0) {
        
        inputGrp <- datalist_Blast_working() %>%
          slice(input$BLAST_rows_selected)
        
        ASVname <- paste(inputGrp$Phylum, inputGrp$Class, inputGrp$Order,
                         inputGrp$Family, inputGrp$Genus, inputGrp$Species, sep = " - ")
        
        plotSingle(datalist_Atlantic %>%
                     makeAppProportion(., Ocean = "Atlantic"), as_vector(select(inputGrp,1)), ASVname)
        
      } else {
        
        ggplot(tibble(text = "ASV not found in Atlantic")) +
          geom_text(aes(x = 1, y = 1, label = text)) +
          theme_void()
        
      }
    }
  })
  
}