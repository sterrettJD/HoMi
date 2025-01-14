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


plot_data <- function(df){
  p <- df %>% 
    pivot_longer(!c(`Percent host`, `Taxonomy method`, project),
                 names_to=c("Taxon"),
                 values_to=c("Abundance")) %>%
    mutate(Abundance=as.numeric(Abundance)) %>%
    ggplot(mapping=aes(x=`Percent host`, y=Abundance, fill=`Taxonomy method`)) +
    geom_boxplot(outliers=F) +
    geom_point(position=position_jitterdodge()) +
    geom_hline(yintercept=1/3) +
    facet_wrap(Taxon ~ project, ncol=2) +
    theme_bw()
  
  return(p)
}


read_file_and_get_level_df <- function(file, level, level.values){
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
  
  # get host percent
  taxon.lvl$`Percent host` <- get_percent_host(rownames(taxon.lvl))
  
  return(taxon.lvl)
}


main <- function(){
  opts <- get_args()
  tax.level <- opts$taxon_level
  split.files <- unlist(str_split(opts$input_files, pattern=","))
  
  if(tax.level=="genus"){
    level.values <- c("Clostridium", "Escherichia", "Faecalibacterium")
  } else if (tax.level=="species"){
    level.values <- c("Clostridium beijerinckii", "Escherichia coli", "Faecalibacterium prausnitzii")
  } else {
    stop("Please pass either `genus` or `species` as the value for --tax_level")
  }
  
  print("Using files:")
  for(file in split.files){
    print(file)
  }
  
  dfs <- lapply(split.files,
         FUN=function(file){
           df <- read_file_and_get_level_df(file, tax.level, level.values)
           
           # Make a taxonomy method column
           taxonomy.method <- infer_kraken_or_metaphlan(file)
           df$`Taxonomy method` <- str_to_title(taxonomy.method)           
           
           # Make each row name specific to the project and taxonomy method 
           # from which it came
           # So we can merge these
           project <- str_split_i(file, "\\.", i=1)
           # hard coding better names
           if(project=="benchmarking_synthetic"){
             project <- "Synthetic communities"
           } else if(project=="benchmarking_synthetic_transcriptomes"){
             project <- "Synthetic transcriptomes"
           }
           df$project <- project
           rownames(df) <- paste0(rownames(df), project, taxonomy.method)
           return(df)
         })
  
  comb.dat <- data.table::rbindlist(dfs)
  
  p <- plot_data(comb.dat)
  ggsave(plot=p, 
         filename=file.path(opts$output_dir, 
                        paste0("combined_taxa_boxplot_",
                               tax.level, ".pdf")),
         height=12, width=12)
}

main()
