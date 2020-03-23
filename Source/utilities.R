### Utlility Scripts for OceanBacteria App ###

plotSingle <- function(datalist, ASVnum, ASVname) {
  # plotSingle() displays the abundance gradient of a single ASV along the complete transect and shows it for different
  # Depth_Grps and Size_Fractions. To choose an ASV enter its #OTU_ID in the ASVnum argument.
  # The #OTU_ID takes a specific ID which stands for a specific nucleotide sequence. 
  
  if (is.character(ASVnum) & sum(as_vector(datalist$Count_Data[,1]) %in% ASVnum) == 1) {
    
    data_selected <- datalist$Count_Data %>%
      filter(.[,1] == ASVnum) %>%
      select_if(is.numeric)
    
  } else {
    
    ASVname <- ""
    datalist$Meta_Data %>%
      mutate(Plot = 0) %>%
      ggplot(.,(aes(x = Latitude, y = Plot))) +
      facet_grid(Depth_Grp~Size_Fraction) +
      labs(title = "") +
      theme(plot.title = element_text(size = 7))
    
    return()
    
  }
  
  data_selected <- data_selected %>%
    select(-Richness) %>%
    t() %>%
    cbind(., datalist$Meta_Data) %>%
    rename("Plot" = ".") 
  
  ggplot(data_selected, aes()) +
    geom_ribbon(aes(x = Latitude, ymin = 0, ymax = Plot), fill = "red") +
    facet_grid(Depth_Grp~Size_Fraction) +
    labs(title = ASVname) +
    theme(plot.title = element_text(size = 7))
  
}

updateTaxonomy <- function(datalist_1, datalist_2) {
  
  Counts_1_df <- datalist_1$Count_Data %>%
    select_if(is.numeric) %>%
    mutate(Counts = rowSums(.)) %>%
    select(Counts, Richness) %>%
    mutate(`#OTU_ID` = as_vector(datalist_1$Count_Data[,1]))
  
  Counts_2_df <- datalist_2$Count_Data %>%
    select_if(is.numeric) %>%
    mutate(Counts = rowSums(.)) %>%
    select(Counts, Richness) %>%
    mutate(`#OTU_ID` = as_vector(datalist_2$Count_Data[,1]))
  
  datalist_Taxonomy <- rbind(select_if(datalist_1$Count_Data, is.character), 
                             select_if(datalist_2$Count_Data, is.character)) %>%
    distinct() %>%
    mutate(Counts_1 = 0, Counts_2 = 0,
           Richness_1 = 0, Richness_2 = 0) %>%
    mutate(Counts_1 = Counts_1_df[match(as_vector(select(.,1)), Counts_1_df$`#OTU_ID`),1] %>%
             as_vector(.) %>%
             remove_NAs(.)) %>%
    mutate(Counts_2 = Counts_2_df[match(as_vector(select(., 1)), Counts_2_df$`#OTU_ID`),1] %>%
             as_vector(.) %>%
             remove_NAs(.)) %>%
    mutate(Richness_1 = Counts_1_df[match(as_vector(select(.,1)), Counts_1_df$`#OTU_ID`),2] %>%
             as_vector(.) %>%
             remove_NAs(.)) %>%
    mutate(Richness_2 = Counts_2_df[match(as_vector(select(.,1)), Counts_2_df$`#OTU_ID`),2] %>%
             as_vector(.) %>%
             remove_NAs(.)) 
  
  return(datalist_Taxonomy)
}

makeAppProportion <- function(datalist, CountsTotal) {
  
  tmp <- datalist$Count_Data %>%
    select_if(is.numeric)
  
  for (i in 1:length(CountsTotal)) {
    tmp <- tmp %>%
      mutate_at(i, function(x) x/CountsTotal[i])
  }
  
  datalist$Count_Data <- datalist$Count_Data %>%
    select_if(is.character) %>%
    cbind(., tmp) %>%
    as_tibble()
  
  return(datalist)
  
}

remove_NAs <- function(Count_Data, replace = 0) {
  
  result <- Count_Data 
  result[is.na(result)] <- replace
  
  return(result)
  
}

combineOceanDatalists <- function(datalist_1, datalist_2) {
  
  datalist_complete <- datalist_1
  
  datalist_complete$Count_Data <- full_join(select(datalist_1$Count_Data, -Richness), 
                                            select(datalist_2$Count_Data, -Richness)) %>%
    remove_NAs()
  
  datalist_complete$Meta_Data <- full_join(mutate(datalist_1$Meta_Data,
                                                  Station = as.character(Station)), 
                                           mutate(datalist_2$Meta_Data,
                                                  Station = as.character(Station)))
  
  return(datalist_complete)
  
}

blastWrapper <- function(blastCommand, blastDB, queryFile) {
  
  colnames <- c("qseqid", "#OTU_ID", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                "send", "evalue", "bitscore")
  
  system2(blastCommand, 
          args = c("-db", blastDB, 
                   "-query", queryFile,
                   "-outfmt", 6, 
                   "-evalue", 1e-6,
                   "-ungapped"), wait = T, stdout = T) %>%
    as_tibble() %>% 
    separate(col = value, 
             into = colnames,
             sep = "\t",
             convert = TRUE) %>%
    select(-c(1,6))
  
}