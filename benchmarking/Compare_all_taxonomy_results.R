if (!require("tidyverse")){
  install.packages("tidyverse", repos="http://cran.us.r-project.org")
  library("tidyverse")
}

if (!require("optparse")){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library("optparse")
}

if (!require("ggplot2")){
  install.packages("ggplot2", repos="http://cran.us.r-project.org")
  library("ggplot2")
}


get_args <- function(){
  option_list <- list( 
    make_option(c("-i", "--input_files"),
                help="The taxonomy abundance tables, comma separated, no spaces"),
    make_option(c("-l", "--taxon_level"),
                help="`genus` or `species`"),
    make_option(c("-r", "--read_lengths"),
                help="A csv file with the average read length per genome.
                There should be two columns, `species`, and `length`."),
    make_option(c("-o", "--output_dir"),
                help="The output directory for plots and results",
                default=NULL)
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  return (opt)
}


infer_kraken_or_metaphlan <- function(filepath){
  if(grepl("kraken", filepath)){
    return("kraken")
  }
  else if (grepl("humann", filepath)){
    return("metaphlan")
  }
  else{
    stop(paste("Could not infer taxonomy method for", filepath))
  }
}

read_taxonomy_data <- function(file, type){
  taxonomy.types <- c("kraken", "metaphlan")
  if((type %in% taxonomy.types)==F){
    stop(paste("Type", type, "is not in options:",
               paste(taxonomy.types, collapse=", ")))
  }
  
  if(type=="kraken"){
    return(read_kraken(file))
  }
  
  if(type=="metaphlan"){
    read_metaphlan(file)
  }
}

read_kraken <- function(file){
  # Need to specify header=TRUE here, which ensures that data.table uses the header.
  # If this is removed, it will not be robust to numeric sample IDs
  # which will trigger data table to not use this row as the column names
  bugslist <- data.table::fread(file, header=TRUE, data.table=F)
  rownames(bugslist) <- bugslist$`tax_id`
  
  # Grab just the feature table
  rel.cols <- colnames(bugslist)[grepl(x=colnames(bugslist), pattern="_frac")]
  feature.table <- bugslist[,rel.cols]
  colnames(feature.table) <- gsub(
    colnames(feature.table), 
    pattern=".bracken_frac", replacement="")
  
  # create tax table
  bugs <- rownames(bugslist)
  tax.columns <- c("domain", "phylum", 
                   "class", "order",
                   "family", "genus",
                   "species")
  taxonomy.table <- bugslist %>%
    select(all_of(tax.columns))
  rownames(taxonomy.table) <- bugs
  
  
  return(list(abundance=feature.table, tax=taxonomy.table))
}


filter_host_from_kraken_df <- function(feature.table, host.id="9606"){
  return(feature.table[rownames(feature.table)!=host.id,])
}


read_metaphlan <- function(file){
  # Need to specify header=TRUE here, which ensures that data.table uses the header.
  # If this is removed, it will not be robust to numeric sample IDs
  # which will trigger data table to not use this row as the column names
  bugslist <- data.table::fread(file, header=TRUE) %>% as.data.frame()
  rownames(bugslist) <- bugslist$`#clade_name`
  bugslist <- bugslist %>% 
    dplyr::select(-c(`#clade_name`,
                     NCBI_tax_id,
                     additional_species))
  
  # just get the lowest tax level. We can collapse again later
  feature.table <- filter_mphlan_by_taxonomy_level(bugslist, level="max")
  
  # replace NA with 0
  feature.table[is.na(feature.table)] <- 0
  
  bugs <- rownames(feature.table)
  taxonomy.table <- metaphlan_names_to_tax_table(bugs)
  rownames(taxonomy.table) <- bugs
  
  return(list(abundance=feature.table, tax=taxonomy.table))
}


# grabs a certain level of taxonomic resolution from the metaphlan output
# defaults to the max level of resolution
filter_mphlan_by_taxonomy_level <- function(df, level="max"){
  tax.names <- rownames(df)
  levels.contained <- stringr::str_count(tax.names, "\\|")
  
  if (level=="max"){
    level <- max(levels.contained)
  }
  
  to.return <- df[levels.contained==level,]
  rownames(to.return) <- rownames(df)[levels.contained==level]
  return (to.return)
}


# Converts metaphlan formatted taxonomy names to a taxonomy table for phyloseq
# based on https://gist.github.com/lwaldron/512d1925a8102e921f05c5b25de7ec94
metaphlan_names_to_tax_table <- function(bugs){
  splitted <- strsplit(bugs, split="|", fixed=T)
  # create empty taxonomy matrix
  taxmat <- matrix(NA, 
                   ncol=max(sapply(splitted, length)), 
                   nrow=length(splitted))
  colnames(taxmat) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")[1:ncol(taxmat)]
  
  # add split taxonomy to the matrix
  for (i in 1:nrow(taxmat)){
    tax.resolution <- length(splitted[[i]])
    taxmat[i, 1:tax.resolution] <- splitted[[i]]
  }
  # remove the p__, f__, etc to indicate level
  taxmat <- gsub("[a-z]__", "", taxmat)
  
  return(as.data.frame(taxmat))
}


calculate_relative_abundance <- function(feature.table){
  column.sums <- colSums(feature.table)
  return(sweep(feature.table, 2, column.sums, FUN = "/"))
}


get_important_abundances <- function(feature.table, tax.table, tax.level, groups){
  features.with.tax <- merge(feature.table, 
                             tax.table, 
                             by="row.names")
  group.summed.table <- features.with.tax %>% 
    dplyr::group_by(!!sym(tax.level)) %>%
    summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), 
              .groups = "drop") %>%
    filter(!!sym(tax.level) %in% groups)
  
  return(group.summed.table)
}


get_percent_host <- function(sample.names){
  
  percents <- sapply(sample.names,
                     FUN=function(x){
                       unlist(strsplit(x, split="_"))[1]}
  )
  return(percents)
}


plot_data <- function(df, level.values){
  
  hline.df <- data.frame(Taxon=c(level.values, 
                                 "Other"),
                         line=c(rep(1/3, length(level.values)), 
                                0)
                         )
  
  p <- df %>% 
    mutate_at(vars(-c(`Percent host`, `Taxonomy method`, project)),
              as.numeric) %>%
    pivot_longer(!c(`Percent host`, `Taxonomy method`, project),
                 names_to=c("Taxon"),
                 values_to=c("Abundance")) %>%
    mutate(Abundance=as.numeric(Abundance)) %>%
    ggplot(mapping=aes(x=`Percent host`, y=Abundance, fill=`Taxonomy method`)) +
    geom_boxplot(outliers=F) +
    geom_point(position=position_jitterdodge()) +
    facet_grid(Taxon ~ project, 
               scales="free") +
    geom_hline(data=hline.df,
               aes(yintercept=line)) +
    theme_bw()
  
  return(p)
}


read_file_and_get_level_df <- function(file, level, level.values, read.lengths, normalize){
  taxonomy.method <- infer_kraken_or_metaphlan(file)
  
  data <- read_taxonomy_data(file, type=taxonomy.method)  
  
  feature.table <- data$abundance
  tax.table <- data$tax
  
  if(taxonomy.method=="kraken"){
    # remove host reads
    feature.table <- filter_host_from_kraken_df(feature.table)
    feature.table <- calculate_relative_abundance(feature.table)
  }
  
  if(taxonomy.method=="metaphlan"){
    feature.table <- calculate_relative_abundance(feature.table)
    tax.table <- tax.table %>%
      mutate_all(function(x){gsub("_", " ", x)})
  }
  
  
  # get feature abundances for the taxa of interest
  # and transpose for plotting
  # and clean up column names
  taxon.lvl <- t(get_important_abundances(feature.table, tax.table, 
                                          level, level.values)) %>%
    as.data.frame()
  colnames(taxon.lvl) <- taxon.lvl[level,]
  taxon.lvl <- taxon.lvl[rownames(taxon.lvl)!=level,]

  # Make a column if it isn't there for each of the taxa that should be detected
  # And make sure each column is numeric
  for(taxon in level.values){
    if((taxon %in% colnames(taxon.lvl))==F){
      taxon.lvl[,taxon] <- 0
    }

    taxon.lvl[,taxon] <- as.numeric(taxon.lvl[,taxon])
  }

  # Normalize taxon by average read length for that taxon
  # Some taxa have shorter reads due to the study they were pulled from
  # So this will be addressed here
  if(normalize==T){
    for(taxon in level.values){
    norm.factor <- read.lengths[read.lengths[level]==taxon, "normalization_factor"]
    taxon.lvl[,taxon] <- taxon.lvl[,taxon]*norm.factor
    }
  }
  
  # get host percent
  taxon.lvl$`Percent host` <- get_percent_host(rownames(taxon.lvl))
  
  return(taxon.lvl)
}


compare_abundances_to_theoretical <- function(df){
  df %>% 
    pivot_longer(!c(`Percent host`, `Taxonomy method`, project),
                 names_to=c("Taxon"),
                 values_to=c("Abundance")) %>%
    mutate(Abundance=as.numeric(Abundance),
           Abundance_diff=Abundance-(1/3)) %>%
    group_by(`Taxonomy method`, project, Taxon) %>% 
    summarise(mean=mean(Abundance_diff),
              sd=sd(Abundance_diff),
              min=min(Abundance_diff),
              max=max(Abundance_diff),
              lower = mean(Abundance_diff) - qt(1 - 0.05/2, (n()-1))*sd(Abundance_diff)/sqrt(n()),
              upper = mean(Abundance_diff) + qt(1 - 0.05/2, (n()-1))*sd(Abundance_diff)/sqrt(n()))
}


get_abundance_not_in_level_values <- function(df, level.values){
  res <- df %>%
    dplyr::select(all_of(level.values)) %>%
    mutate_all(as.numeric) %>%
    rowSums()
  return(1-res)
}

summarize_abundance_not_in_level_values <- function(df){
  res <- df %>%
    # The values should be 0 here, which means that abund = the diff
    mutate(Abundance=as.numeric(`Other`),
           Abundance_diff=Abundance-(0)) %>%
    group_by(`Taxonomy method`, project) %>% 
    summarise(mean=mean(Abundance_diff),
              sd=sd(Abundance_diff),
              min=min(Abundance_diff),
              max=max(Abundance_diff),
              lower = mean(Abundance_diff) - qt(1 - 0.05/2, (n()-1))*sd(Abundance_diff)/sqrt(n()),
              upper = mean(Abundance_diff) + qt(1 - 0.05/2, (n()-1))*sd(Abundance_diff)/sqrt(n()))
    res$Taxon <- "Other"
    return(res)
}


plot_abundance_not_in_level_values <- function(df){
  p <- df %>%
    ggplot(mapping=aes(x=`Percent host`, y=`Other`, 
                       fill=`Taxonomy method`)) +
    geom_boxplot(outliers=F) +
    geom_point(position=position_jitterdodge()) + 
    facet_grid( ~ project) +
    theme_bw()
    
  return(p)
}

main <- function(){
  opts <- get_args()
  tax.level <- opts$taxon_level
  split.files <- unlist(str_split(opts$input_files, pattern=","))
  
  if(tax.level=="genus"){
    level.values <- c("Clostridium", "Escherichia", "Bacteroides")
  } else if (tax.level=="species"){
    level.values <- c("Clostridium beijerinckii", "Escherichia coli", "Bacteroides fragilis")
  } else {
    stop("Please pass either `genus` or `species` as the value for --tax_level")
  }
  

  print(paste("Reading read length file:", opts$read_lengths))
  read.lengths <- data.table::fread(opts$read_lengths, data.table=F)
  if(!setequal(colnames(read.lengths), c("species", "length"))){
    stop("Column names for read_lengths should be `species` and `length`")
  }
  read.lengths <- filter(read.lengths, species!="human")
  read.lengths$genus <- str_split_i(read.lengths$species, " ", i=1)
  read.lengths.total <- sum(read.lengths$length)
  average.read.length <- read.lengths.total/nrow(read.lengths)
  read.lengths$normalization_factor <- average.read.length/read.lengths$length


  print("Using files:")
  for(file in split.files){
    print(file)
  }
  
  dfs <- lapply(split.files,
         FUN=function(file){
          
           project <- str_split_i(file, "\\.", i=1)
           
           # hard coding better project names
           if(project %in% c("dna_benchmarking_synthetic", 
                             "rna_benchmarking_synthetic")){
             project <- "Synthetic communities"
             norm.by.read.length <- FALSE
           } else if(project %in% c("dna_benchmarking_synthetic_transcriptomes",
                                    "rna_benchmarking_synthetic_transcriptomes")){
             project <- "Synthetic transcriptomes"
             norm.by.read.length <- FALSE
           } else if(project %in% c("dna_semi", "rna_semi")){
            project <- "Semisynthetic transcriptomes"
            norm.by.read.length <- TRUE
           } else {
            stop(paste("Project name", project, "isn't recognized."))
           }
           # get the data
           df <- read_file_and_get_level_df(file, tax.level, level.values, 
                                            read.lengths, normalize=norm.by.read.length)
           
           # Make a taxonomy method column
           taxonomy.method <- infer_kraken_or_metaphlan(file)
           df$`Taxonomy method` <- str_to_title(taxonomy.method)           
           
           # Make each row name specific to the project and taxonomy method 
           # from which it came
           # So we can merge these
           df$project <- project
           rownames(df) <- paste0(rownames(df), project, taxonomy.method)
           return(df)
         })

  for(df in dfs){
    print(head(df))
  }
  comb.dat <- data.table::rbindlist(dfs)
  compare.df <- compare_abundances_to_theoretical(comb.dat)  
  
  # Get the "Other" abundances here
  comb.dat$`Other` <- get_abundance_not_in_level_values(comb.dat, 
                                                          level.values)
  p <- plot_abundance_not_in_level_values(comb.dat)
  ggsave(plot=p, 
         filename=file.path(opts$output_dir, 
                            paste0("unassigned_boxplot_",
                                   tax.level, ".pdf")))
  # Add the "other" values
  compare.df <- rbind(summarize_abundance_not_in_level_values(comb.dat), 
                      compare.df)
  compare.df %>%
    mutate(`mean [95% CI]` = paste0(round(mean, digits=3), 
                                    " [", round(lower, digits=3),
                                    " - ", round(upper, digits=3), "]"),
           min = round(min, 3),
           max = round(max, 3)) %>%
    select(c(`Taxonomy method`, Taxon, project, `mean [95% CI]`, min, max)) %>%
    pivot_wider(names_from = `Taxonomy method`, 
                values_from = c(`mean [95% CI]`, max, min),
                names_glue = "{`Taxonomy method`} {.value}"
    ) %>% 
    write.csv(file=file.path(opts$output_dir, 
                              paste0("summarized_benchmark_",
                              tax.level, ".csv")))
  
  p <- plot_data(comb.dat, level.values)
  ggsave(plot=p, 
         filename=file.path(opts$output_dir, 
                            paste0("combined_taxa_boxplot_",
                                   tax.level, ".pdf")),
         height=12, width=12)
}

main()
