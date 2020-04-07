### Serverside output rendering ###

server = function(input, output, session) {
  
  #### Prokaryotes ####
  observeEvent(input$invBlastGo, {
    if (grepl("^[ATCG]+$", input$hashInput)) {
      showNotification("Enter OTU_IDs. Sequences not allowed.")
      output$hashOutput <- renderText({""})
    } else {
      output$hashOutput <- renderText({
        system(paste("grep -A 1", input$hashInput, "Data/BlastDB/OceanBacteriaSequences.fasta"), intern=T)[2]
      })
    }
  })
  
  observeEvent(input$ASV_rows_selected, {
    
    if (is.numeric(input$ASV_rows_selected)) {
      output$MapOutput <- renderLeaflet({
        
        tmp <- provinces
        tmp@data <- Meta_Combined_working() %>%
          filter(Size_Fraction == ifelse(input$MapSizeFraction == "0.22-3µm (FL)", "0.22",
                                  ifelse(input$MapSizeFraction == "3-8µm (SPA)", "3",
                                  ifelse(input$MapSizeFraction == ">8µm (LPA)", "8", "")))) %>%
          filter(Depth_Grp == input$MapDepth) %>%
          mutate(Province = ifelse(Province == "NAST", "NASE", 
                            ifelse(Province == "NTPG", "NPTG",
                            ifelse(Province == "NPST", "NPSW",
                            ifelse(Province == "PSAG", "PSAE", Province))))) %>%
          group_by(Province) %>%
          summarize(Abundance = mean(Abundance)) %>%
          rename("ProvCode" = "Province") %>%
          left_join(provinces@data, ., by = "ProvCode")
        
        leaflet(tmp) %>%
            addTiles() %>%
            #addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
            addPolygons(weight = 1, color = "#444444",smoothFactor = 0.5,
                      opacity = 1, fillOpacity = 0.5, 
                      fillColor = ~colorNumeric("YlOrRd", c(input$ColorRange[1], max(input$ColorRange[2], max(Abundance, na.rm = T))))(Abundance))  %>%
            fitBounds(lng1 = -70, lat1 = -77, lng2 = 70, lat2 = 80) %>%
          addLabelOnlyMarkers(data = province_center, lng = ~long, lat = ~lat, label = ~ProvCode,
                              labelOptions = labelOptions(noHide = T, direction = 'top', textOnly = T))
      
      })
    }
  })
  
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
                                       `OTU_ID` == {datalist_Blast_working() %>% 
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
    
    query <- input$seqInputBlast %>%
      gsub("[\n\r ]", "", .)
    
    if (is.character(query) & nchar(query) > 4 &
        grepl("^[ATCGYRWSKMDVHBXN]+$", query)) {
      
      write_file(paste(">Query\n", input$seqInputBlast, sep = ""), queryFile)
      
      blastTable <- blastWrapper(blastCommand, blastDB, queryFile) %>%
        filter(pident > 0.95)
      
      output <- right_join(datalist_Taxonomy, blastTable, by = "OTU_ID") %>%
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
    selected = "All",
    options = list(placeholder = 'select'),
    server = TRUE
  )
  
  observeEvent(input$ASV_rows_selected | input$BLAST_rows_selected, {
    updateSelectInput(
      session = session,
      inputId = "Xparam",
      choices = names(Meta_Combined_working()),
      selected = "Pot_Temperature"
      )
    updateSelectInput(
      session = session,
      inputId = "Yparam",
      choices = names(Meta_Combined_working()),
      selected = "Abundance"
    )
  })

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
                    mutate(Counts_1 = round(Counts_1 / sum(sum(select_if(datalist_Atlantic$Count_Data, is.numeric))), digits = 3) * 100,
                           Counts_2 = round(Counts_2 / sum(sum(select_if(datalist_Pacific$Count_Data, is.numeric))), digits = 3) * 100) %>%
                    rename("Proportion Atlantic" = "Counts_1",
                           "Proportion Pacific" = "Counts_2",
                           "Richness Atlantic" = "Richness_1",
                           "Richness Pacific" = "Richness_2"),
                  selection = "single")
  )
  
  output$BLAST <- DT::renderDataTable(
    DT::datatable(datalist_Blast_working() %>%
                    mutate(Counts_1 = round(Counts_1 / sum(sum(select_if(datalist_Atlantic$Count_Data, is.numeric))), digits = 3) * 100,
                           Counts_2 = round(Counts_2 / sum(sum(select_if(datalist_Pacific$Count_Data, is.numeric))), digits = 3) * 100) %>%
                    rename("Proportion Atlantic" = "Counts_1",
                           "Proportion Pacific" = "Counts_2",
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
                     makeAppProportion(., apply(select_if(datalist_Pacific$Count_Data, is.numeric), 2, sum)), inputGrp, ASVname)
        
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
                     makeAppProportion(., apply(select_if(datalist_Atlantic$Count_Data, is.numeric), 2, sum)), inputGrp, ASVname)
        
        
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
                     makeAppProportion(., apply(select_if(datalist_Pacific$Count_Data, is.numeric), 2, sum)), as_vector(select(inputGrp,1)), ASVname)
        
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
                     makeAppProportion(., apply(select_if(datalist_Atlantic$Count_Data, is.numeric), 2, sum)), as_vector(select(inputGrp,1)), ASVname)
        
      } else {
        
        ggplot(tibble(text = "ASV not found in Atlantic")) +
          geom_text(aes(x = 1, y = 1, label = text)) +
          theme_void()
        
      }
    }
  })
  
  #### Chloroplasts ####
  
  observeEvent(input$ChloroInvBlastGo, {
    if (grepl("^[ATCG]+$", input$ChloroHashInput)) {
      showNotification("Enter OTU_IDs. Sequences not allowed.")
      output$ChloroHashOutput <- renderText({""})
    } else {
      output$ChloroHashOutput <- renderText({
        system(paste("grep -A 1", input$ChloroHashInput, "Data/BlastDB/OceanChloroplastSequences.fasta"), intern=T)[2]
      })
    }
  })
  
  observeEvent(input$ChloroASV_rows_selected, {
    
    if (is.numeric(input$ChloroASV_rows_selected)) {
      output$ChloroMapOutput <- renderLeaflet({
        
        tmp <- provinces
        tmp@data <- Chloro_Meta_Combined_working() %>%
          filter(Size_Fraction == ifelse(input$ChloroMapSizeFraction == "0.22-3µm (FL)", "0.22",
                                  ifelse(input$ChloroMapSizeFraction == "3-8µm (SPA)", "3",
                                  ifelse(input$ChloroMapSizeFraction == ">8µm (LPA)", "8", "")))) %>%
          filter(Depth_Grp == input$ChloroMapDepth) %>%
          mutate(Province = ifelse(Province == "NAST", "NASE", 
                            ifelse(Province == "NTPG", "NPTG",
                            ifelse(Province == "NPST", "NPSW",
                            ifelse(Province == "PSAG", "PSAE", Province))))) %>%
          group_by(Province) %>%
          summarize(Abundance = mean(Abundance)) %>%
          rename("ProvCode" = "Province") %>%
          left_join(provinces@data, ., by = "ProvCode")
        
        leaflet(tmp) %>%
          addTiles() %>%
          addPolygons(weight = 1, color = "#444444",smoothFactor = 0.5,
                      opacity = 1, fillOpacity = 0.5, 
                      fillColor = ~colorNumeric("YlOrRd", c(input$ChloroColorRange[1], max(input$ChloroColorRange[2], max(Abundance, na.rm = T))))(Abundance))  %>%
          fitBounds(lng1 = -70, lat1 = -77, lng2 = 70, lat2 = 80) %>%
          addLabelOnlyMarkers(data = province_center, lng = ~long, lat = ~lat, label = ~ProvCode,
                              labelOptions = labelOptions(noHide = T, direction = 'top', textOnly = T))
        
      })
    }
  })
  
  output$ChloroPlotCorAnalysis <- renderPlot({
    
    if (input$ChloroColparam == "Empty" | is_empty(input$ChloroColparam) | input$ChloroColparam == "") {
      colorValue <- "red"
    } else {
      colorValue <- input$ChloroColparam
    }
    
    if (input$ChloroSizeparam == "Empty" | is_empty(input$ChloroSizeparam) | input$ChloroSizeparam == "") {
      sizeValue <- 10
    } else {
      sizeValue <- input$ChloroSizeparam
    }
    
    Chloro_Meta_Combined_working() %>%
      ggplot(., aes(x = !!input$ChloroXparam, y = !!input$ChloroYparam, col = !!input$ChloroColparam, size = !!input$ChloroSizeparam)) +
      geom_point()
    
  })
  
  Chloro_Meta_Combined_working <- reactive({
    
    if (is.numeric(input$ChloroASV_rows_selected)) {
      
      Chloro_datalist_Combined <- combineOceanDatalists(
        Chloro_datalist_Atlantic_working()  %>%
          makeAppProportion(., apply(select_if(Chloro_datalist_Atlantic$Count_Data, is.numeric), 2, sum)), 
        Chloro_datalist_Pacific_working()  %>%
          makeAppProportion(., apply(select_if(Chloro_datalist_Pacific$Count_Data, is.numeric), 2, sum))
      )
      Chloro_Meta_Combined_working <- Chloro_datalist_Combined$Meta_Data %>%
        mutate(Abundance = Chloro_datalist_Combined$Count_Data %>%
                 slice(input$ChloroASV_rows_selected) %>%
                 select_if(is.numeric) %>%
                 as.numeric()
        )
      
    } else if (is.numeric(input$ChloroBLAST_rows_selected)) {
      
      Chloro_datalist_Combined <- combineOceanDatalists(
        Chloro_datalist_Atlantic %>%
          makeAppProportion(., apply(select_if(Chloro_datalist_Atlantic$Count_Data, is.numeric), 2, sum)), 
        Chloro_datalist_Pacific %>%
          makeAppProportion(., apply(select_if(Chloro_datalist_Pacific$Count_Data, is.numeric), 2, sum))
      )
      Chloro_Meta_Combined_working <- Chloro_datalist_Combined$Meta_Data %>%
        mutate(Abundance = Chloro_datalist_Combined$Count_Data %>%
                 filter(
                   `OTU_ID` == {Chloro_datalist_Blast_working() %>% 
                       slice(input$ChloroBLAST_rows_selected) %>%
                       select(1) %>%
                       as_vector()
                   }) %>%
                 select_if(is.numeric) %>%
                 slice(1) %>%
                 as.numeric()
        )
      
    } else {
      Chloro_datalist_Combined$Meta_Data
    }
  })
  
  Chloro_datalist_Atlantic_working <- reactive({
    
    tmp <- Chloro_datalist_Atlantic
    
    if (input$ChloroAbundanceFilter) {
      tmp <- tmp %>%
        filterMiliciStyle(.)
    }
    
    if (input$ChloroSupergroup == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        group_by(Supergroup) %>%
        summarize_if(is.numeric, sum)
      
    } else if (input$ChloroPhylum == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        group_by(Phylum) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroClass == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        group_by(Class) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroSubclass == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        group_by(Subclass) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroOrder == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        group_by(Order) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroSuborder == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        group_by(Suborder) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroFamily == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        filter(Suborder == input$ChloroSuborder) %>%
        group_by(Family) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroGenus == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        filter(Suborder == input$ChloroSuborder) %>%
        filter(Family == input$ChloroFamily) %>%
        group_by(Genus) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroSpecies == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        filter(Suborder == input$ChloroSuborder) %>%
        filter(Family == input$ChloroFamily) %>%
        filter(Genus == input$ChloroGenus) %>%
        group_by(Species) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        filter(Suborder == input$ChloroSuborder) %>%
        filter(Family == input$ChloroFamily) %>%
        filter(Genus == input$ChloroGenus) %>%
        filter(Species == input$ChloroSpecies)
      
    }
    tmp
  })
  
  Chloro_datalist_Pacific_working <- reactive({
    
    tmp <- Chloro_datalist_Pacific
    
    if (input$ChloroAbundanceFilter) {
      tmp <- tmp %>%
        filterMiliciStyle(.)
    }
    
    if (input$ChloroSupergroup == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        group_by(Supergroup) %>%
        summarize_if(is.numeric, sum)
      
    } else if (input$ChloroPhylum == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        group_by(Phylum) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroClass == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        group_by(Class) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroSubclass == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        group_by(Subclass) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroOrder == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        group_by(Order) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroSuborder == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        group_by(Suborder) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroFamily == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        filter(Suborder == input$ChloroSuborder) %>%
        group_by(Family) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroGenus == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        filter(Suborder == input$ChloroSuborder) %>%
        filter(Family == input$ChloroFamily) %>%
        group_by(Genus) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else if (input$ChloroSpecies == "All") {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        filter(Suborder == input$ChloroSuborder) %>%
        filter(Family == input$ChloroFamily) %>%
        filter(Genus == input$ChloroGenus) %>%
        group_by(Species) %>%
        select_if(is.numeric) %>%
        summarize_all(sum)
      
    } else {
      
      tmp$Count_Data <- tmp$Count_Data %>%
        filter(Supergroup == input$ChloroSupergroup) %>%
        filter(Phylum == input$ChloroPhylum) %>%
        filter(Class == input$ChloroClass) %>%
        filter(Subclass == input$ChloroSubclass) %>%
        filter(Order == input$ChloroOrder) %>%
        filter(Suborder == input$ChloroSuborder) %>%
        filter(Family == input$ChloroFamily) %>%
        filter(Genus == input$ChloroGenus) %>%
        filter(Species == input$ChloroSpecies)
      
    }
    tmp
  })
  
  Chloro_datalist_Taxonomy_working <- reactive({ 
    
    updateTaxonomy(Chloro_datalist_Atlantic_working(), Chloro_datalist_Pacific_working())
    
  })
  
  Chloro_datalist_Blast_working <- eventReactive(input$ChloroBlastGo, {
    
    query <- input$ChloroSeqInputBlast %>%
      gsub("[\n\r ]", "", .)
    
    if (is.character(query) & nchar(query) > 4 &
        grepl("^[ATCGYRWSKMDVHBXN]+$", query)) {
      
      write_file(paste(">Query\n", input$ChloroSeqInputBlast, sep = ""), queryFile)
      
      blastTable <- blastWrapper(blastCommand, ChloroBlastDB, queryFile) %>%
        filter(pident > 0.95)
      
      output <- right_join(Chloro_datalist_Taxonomy, blastTable, by = "OTU_ID") %>%
        drop_na()
      
      output
      
    } else {
      
      showNotification("No sequence input for blast or sequence is too short (min length = 4).")
      tibble()
      
    }
  })
  
  updateSelectizeInput(
    session = session,
    inputId = "ChloroSupergroup",
    choices = c("All", Chloro_datalist_Taxonomy %>%
                  select(Supergroup) %>%
                  as_vector() %>%
                  unique() %>%
                  sort()
    ),
    selected = "All",
    options = list(placeholder = 'select'),
    server = TRUE
  )
  
  observeEvent(input$ChloroASV_rows_selected | input$ChloroBLAST_rows_selected, {
    updateSelectInput(
      session = session,
      inputId = "ChloroXparam",
      choices = names(Chloro_Meta_Combined_working()),
      selected = "Pot_Temperature"
    )
    updateSelectInput(
      session = session,
      inputId = "ChloroYparam",
      choices = names(Chloro_Meta_Combined_working()),
      selected = "Abundance"
    )
  })
  
  observeEvent(input$ChloroSupergroup, {
    choice_Phylum <- Chloro_datalist_Taxonomy %>%
      filter(if (input$ChloroSupergroup != "All") {.$Supergroup == input$ChloroSupergroup} else {TRUE}) %>%
      select(Phylum) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "ChloroPhylum",
      choices = c("All", choice_Phylum))
    
  })
  
  observeEvent(input$ChloroPhylum, {
    choice_Class <- Chloro_datalist_Taxonomy %>%
      filter(if (input$ChloroSupergroup != "All") {.$Supergroup == input$ChloroSupergroup} else {TRUE}) %>%
      filter(if (input$ChloroPhylum != "All") {.$Phylum == input$ChloroPhylum} else {TRUE}) %>%
      select(Class) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "ChloroClass",
      choices = c("All", choice_Class)
    )
  })
  
  observeEvent(input$ChloroClass, {
    choice_Subclass <- Chloro_datalist_Taxonomy %>%
      filter(if (input$ChloroSupergroup != "All") {.$Supergroup == input$ChloroSupergroup} else {TRUE}) %>%
      filter(if (input$ChloroPhylum != "All") {.$Phylum == input$ChloroPhylum} else {TRUE}) %>%
      filter(if (input$ChloroClass != "All") {.$Class == input$ChloroClass} else {TRUE}) %>%
      select(Subclass) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "ChloroSubclass",
      choices = c("All", choice_Subclass)
    )
  })
  
  observeEvent(input$ChloroSubclass, {
    choice_Order <- Chloro_datalist_Taxonomy %>%
      filter(if (input$ChloroSupergroup != "All") {.$Supergroup == input$ChloroSupergroup} else {TRUE}) %>%
      filter(if (input$ChloroPhylum != "All") {.$Phylum == input$ChloroPhylum} else {TRUE}) %>%
      filter(if (input$ChloroClass != "All") {.$Class == input$ChloroClass} else {TRUE}) %>%
      filter(if (input$ChloroSubclass != "All") {.$Subclass == input$ChloroSubclass} else {TRUE}) %>%
      select(Order) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "ChloroOrder",
      choices = c("All", choice_Order)
    )
  })
  
  observeEvent(input$ChloroOrder, {
    choice_Suborder <- Chloro_datalist_Taxonomy %>%
      filter(if (input$ChloroSupergroup != "All") {.$Supergroup == input$ChloroSupergroup} else {TRUE}) %>%
      filter(if (input$ChloroPhylum != "All") {.$Phylum == input$ChloroPhylum} else {TRUE}) %>%
      filter(if (input$ChloroClass != "All") {.$Class == input$ChloroClass} else {TRUE}) %>%
      filter(if (input$ChloroSubclass != "All") {.$Subclass == input$ChloroSubclass} else {TRUE}) %>%
      filter(if (input$ChloroOrder != "All") {.$Order == input$ChloroOrder} else {TRUE}) %>%
      select(Suborder) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "ChloroSuborder",
      choices = c("All", choice_Suborder)
    )
  })
  
  observeEvent(input$ChloroSuborder, {
    choice_Family <- Chloro_datalist_Taxonomy %>%
      filter(if (input$ChloroSupergroup != "All") {.$Supergroup == input$ChloroSupergroup} else {TRUE}) %>%
      filter(if (input$ChloroPhylum != "All") {.$Phylum == input$ChloroPhylum} else {TRUE}) %>%
      filter(if (input$ChloroClass != "All") {.$Class == input$ChloroClass} else {TRUE}) %>%
      filter(if (input$ChloroSubclass != "All") {.$Subclass == input$ChloroSubclass} else {TRUE}) %>%
      filter(if (input$ChloroOrder != "All") {.$Order == input$ChloroOrder} else {TRUE}) %>%
      filter(if (input$ChloroSuborder != "All") {.$Suborder == input$ChloroSuborder} else {TRUE}) %>%
      select(Family) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "ChloroFamily",
      choices = c("All", choice_Family)
    )
  })
  
  observeEvent(input$ChloroFamily, {
    choice_Genus <- Chloro_datalist_Taxonomy %>%
      filter(if (input$ChloroSupergroup != "All") {.$Supergroup == input$ChloroSupergroup} else {TRUE}) %>%
      filter(if (input$ChloroPhylum != "All") {.$Phylum == input$ChloroPhylum} else {TRUE}) %>%
      filter(if (input$ChloroClass != "All") {.$Class == input$ChloroClass} else {TRUE}) %>%
      filter(if (input$ChloroSubclass != "All") {.$Subclass == input$ChloroSubclass} else {TRUE}) %>%
      filter(if (input$ChloroOrder != "All") {.$Order == input$ChloroOrder} else {TRUE}) %>%
      filter(if (input$ChloroSuborder != "All") {.$Suborder == input$ChloroSuborder} else {TRUE}) %>%
      filter(if (input$ChloroFamily != "All") {.$Family == input$ChloroFamily} else {TRUE}) %>%
      select(Genus) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "ChloroGenus",
      choices = c("All", choice_Genus)
    )
  })
  
  observeEvent(input$ChloroGenus, {
    choice_Species <- Chloro_datalist_Taxonomy %>%
      filter(if (input$ChloroSupergroup != "All") {.$Supergroup == input$ChloroSupergroup} else {TRUE}) %>%
      filter(if (input$ChloroPhylum != "All") {.$Phylum == input$ChloroPhylum} else {TRUE}) %>%
      filter(if (input$ChloroClass != "All") {.$Class == input$ChloroClass} else {TRUE}) %>%
      filter(if (input$ChloroSubclass != "All") {.$Subclass == input$ChloroSubclass} else {TRUE}) %>%
      filter(if (input$ChloroOrder != "All") {.$Order == input$ChloroOrder} else {TRUE}) %>%
      filter(if (input$ChloroSuborder != "All") {.$Suborder == input$ChloroSuborder} else {TRUE}) %>%
      filter(if (input$ChloroFamily != "All") {.$Family == input$ChloroFamily} else {TRUE}) %>%
      filter(if (input$ChloroGenus != "All") {.$Genus == input$ChloroGenus} else {TRUE}) %>%
      select(Species) %>%
      as_vector() %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session = session,
      inputId = "ChloroSpecies",
      choices = c("All", choice_Species)
    )
  })
  
  output$ChloroASV <- DT::renderDataTable(
    DT::datatable(Chloro_datalist_Taxonomy_working() %>%
                    mutate(Counts_1 = round(Counts_1 / sum(sum(select_if(Chloro_datalist_Atlantic$Count_Data, is.numeric))), digits = 3) * 100,
                           Counts_2 = round(Counts_2 / sum(sum(select_if(Chloro_datalist_Pacific$Count_Data, is.numeric))), digits = 3) * 100) %>%
                    rename("Proportion Atlantic" = "Counts_1",
                           "Proportion Pacific" = "Counts_2",
                           "Richness Atlantic" = "Richness_1",
                           "Richness Pacific" = "Richness_2"),
                  selection = "single")
  )
  
  output$ChloroBLAST <- DT::renderDataTable(
    DT::datatable(Chloro_datalist_Blast_working() %>%
                    mutate(Counts_1 = round(Counts_1 / sum(sum(select_if(Chloro_datalist_Atlantic$Count_Data, is.numeric))), digits = 3) * 100,
                           Counts_2 = round(Counts_2 / sum(sum(select_if(Chloro_datalist_Pacific$Count_Data, is.numeric))), digits = 3) * 100) %>%
                    rename("Proportion Atlantic" = "Counts_1",
                           "Proportion Pacific" = "Counts_2",
                           "Richness Atlantic" = "Richness_1",
                           "Richness Pacific" = "Richness_2") %>%
                    makeAppProportion(., ),
                  selection = "single")
  )
  
  output$ChloroPacificBrowse <- renderPlot({
    
    if (is.numeric(input$ChloroASV_rows_selected)) {
      
      if (select(slice(Chloro_datalist_Taxonomy_working(), input$ChloroASV_rows_selected), Counts_2) > 0) {
        
        inputGrp <- Chloro_datalist_Taxonomy_working() %>%
          slice(input$ChloroASV_rows_selected) %>%
          select(1) %>%
          as_vector()
        
        ASVname <- paste(ifelse(input$ChloroSupergroup == "All", inputGrp, input$ChloroSupergroup),
                         ifelse(input$ChloroPhylum == "All", inputGrp, input$ChloroPhylum),
                         ifelse(input$ChloroClass == "All", inputGrp, input$ChloroClass),
                         ifelse(input$ChloroSubclass == "All", inputGrp, input$ChloroSubclass),
                         ifelse(input$ChloroOrder == "All", inputGrp, input$ChloroOrder),
                         ifelse(input$ChloroSuborder == "All", inputGrp, input$ChloroSuborder),
                         ifelse(input$ChloroFamily == "All", inputGrp, input$ChloroFamily),
                         ifelse(input$ChloroGenus == "All", inputGrp, input$ChloroGenus),
                         ifelse(input$ChloroSpecies == "All", inputGrp, input$ChloroSpecies), 
                         sep = " - ")
        
        plotSingle(Chloro_datalist_Pacific_working() %>%
                     makeAppProportion(., apply(select_if(Chloro_datalist_Pacific$Count_Data, is.numeric), 2, sum)), inputGrp, ASVname)
        
      } else {
        
        ggplot(tibble(text = "ASV not found in Pacific.")) +
          geom_text(aes(x = 1, y = 1, label = text)) +
          theme_void()
        
      }
    }
  })
  
  output$ChloroAtlanticBrowse <- renderPlot({
    
    if (is.numeric(input$ChloroASV_rows_selected)) {
      
      if (select(slice(Chloro_datalist_Taxonomy_working(), input$ChloroASV_rows_selected), Counts_1) > 0) {
        
        inputGrp <- Chloro_datalist_Taxonomy_working() %>%
          slice(input$ChloroASV_rows_selected) %>%
          select(1) %>%
          as_vector()
        
        ASVname <- paste(ifelse(input$ChloroSupergroup == "All", inputGrp, input$ChloroSupergroup),
                         ifelse(input$ChloroPhylum == "All", inputGrp, input$ChloroPhylum),
                         ifelse(input$ChloroClass == "All", inputGrp, input$ChloroClass),
                         ifelse(input$ChloroSubclass == "All", inputGrp, input$ChloroSubclass),
                         ifelse(input$ChloroOrder == "All", inputGrp, input$ChloroOrder),
                         ifelse(input$ChloroSuborder == "All", inputGrp, input$ChloroSuborder),
                         ifelse(input$ChloroFamily == "All", inputGrp, input$ChloroFamily),
                         ifelse(input$ChloroGenus == "All", inputGrp, input$ChloroGenus),
                         ifelse(input$ChloroSpecies == "All", inputGrp, input$ChloroSpecies), 
                         sep = " - ")
        
        plotSingle(Chloro_datalist_Atlantic_working() %>%
                     makeAppProportion(., apply(select_if(Chloro_datalist_Atlantic$Count_Data, is.numeric), 2, sum)), inputGrp, ASVname)
        
      } else {
        
        ggplot(tibble(text = "ASV not found in Atlantic")) +
          geom_text(aes(x = 1, y = 1, label = text)) +
          theme_void()
        
      }
    }
  })
  
  output$ChloroPacificBlast <- renderPlot({
    
    if (is.numeric(input$ChloroBLAST_rows_selected)) {
      
      if (select(slice(Chloro_datalist_Blast_working(), input$ChloroBLAST_rows_selected), Counts_2) > 0) {
        
        inputGrp <- Chloro_datalist_Blast_working() %>%
          slice(input$ChloroBLAST_rows_selected)
        
        ASVname <- paste(inputGrp$Supergroup, inputGrp$Phylum, inputGrp$Class, inputGrp$Subclass,
                         inputGrp$Order, inputGrp$Suborder, inputGrp$Family, inputGrp$Genus,
                         inputGrp$Species, sep = " - ")
        
        plotSingle(Chloro_datalist_Pacific %>%
                     makeAppProportion(., apply(select_if(Chloro_datalist_Pacific$Count_Data, is.numeric), 2, sum)), as_vector(select(inputGrp,1)), ASVname)
        
      } else {
        
        ggplot(tibble(text = "ASV not found in Pacific.")) +
          geom_text(aes(x = 1, y = 1, label = text)) +
          theme_void()
        
      }
    }
  })
  
  output$ChloroAtlanticBlast <- renderPlot({
    
    if (is.numeric(input$ChloroBLAST_rows_selected)) {
      
      if (select(slice(Chloro_datalist_Blast_working(), input$ChloroBLAST_rows_selected), Counts_1) > 0) {
        
        inputGrp <- Chloro_datalist_Blast_working() %>%
          slice(input$ChloroBLAST_rows_selected)
        
        ASVname <- paste(inputGrp$Supergroup, inputGrp$Phylum, inputGrp$Class, inputGrp$Subclass,
                         inputGrp$Order, inputGrp$Suborder, inputGrp$Family, inputGrp$Genus,
                         inputGrp$Species, sep = " - ")
        
        plotSingle(Chloro_datalist_Atlantic %>%
                     makeAppProportion(., apply(select_if(Chloro_datalist_Atlantic$Count_Data, is.numeric), 2, sum)), as_vector(select(inputGrp,1)), ASVname)
        
      } else {
        
        ggplot(tibble(text = "ASV not found in Atlantic")) +
          geom_text(aes(x = 1, y = 1, label = text)) +
          theme_void()
        
      }
    }
  })
  
}