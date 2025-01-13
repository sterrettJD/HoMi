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
    make_option(c("-i", "--input_file"),
                help="The taxonomy abundance table"),
    make_option(c("-t", "--taxonomy_method"),
                help="Is this a kraken or metaphlan file?
                      Options include `kraken` and `metaphlan`"),
    make_option(c("--sample_names_are_percent"),
                help="Sample names are already percent host
                      and don't need to be parsed into percent host.
                      Pass this for Pereira",
                default=FALSE),
    make_option(c("-o", "--output_dir"),
                help="The output directory for plots and results",
                default=NULL)
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  return (opt)
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
  # Try to convert directly to percent
  # This works for the pereira dataset
  direct.to.percents <- as.numeric(sample.names)
  # Check if all sample names can be directly converted
  # The synthetic samples will not be able to do this 
  # and should return NA
  # if it works, return. Otherwise, parse sample names
  if(sum(is.na(direct.to.percents))==0){
    return(direct.to.percents)
  }
  
  # parse the names according to the schema for synthetic samples
  percents <- sapply(sample.names,
         FUN=function(x){
           unlist(strsplit(x, split="_"))[1]}
         )
  return(percents)
}


plot_data <- function(df){
  p <- df %>% 
    pivot_longer(!`Percent host`,
                 names_to=c("Taxon"),
                 values_to=c("Abundance")) %>%
    mutate(Abundance=as.numeric(Abundance)) %>%
    ggplot(mapping=aes(x=`Percent host`, y=Abundance)) +
    geom_boxplot(outliers=F) +
    geom_jitter() +
    facet_wrap(vars(Taxon), ncol=1) +
    theme_bw()
  
  return(p)
}


main <- function(){
  opts <- get_args()
  data <- read_taxonomy_data(opts$input_file, type=opts$taxonomy_method)  
  
  feature.table <- data$abundance
  tax.table <- data$tax
  
  genera <- c("Clostridium", "Escherichia", "Faecalibacterium")
  species <- c("Clostridium beijerinckii", "Escherichia coli", "Faecalibacterium prausnitzii")
  
  if(opts$taxonomy_method=="kraken"){
    # remove host reads
    feature.table <- filter_host_from_kraken_df(feature.table)
    feature.table <- calculate_relative_abundance(feature.table)
  }
  
  if(opts$taxonomy_method=="metaphlan"){
    feature.table <- calculate_relative_abundance(feature.table)
    tax.table <- tax.table %>%
      mutate_all(function(x){gsub("_", " ", x)})
  }
  
  # get feature abundances for the taxa of interest
  # and transpose for plotting
  # and clean up column names
  genus.lvl <- t(get_important_abundances(feature.table, tax.table, 
                                          "genus", genera)) %>%
    as.data.frame()
  colnames(genus.lvl) <- genus.lvl["genus",]
  genus.lvl <- genus.lvl[rownames(genus.lvl)!="genus",]
  
  species.lvl <- t(get_important_abundances(feature.table, tax.table, 
                                            "species", species)) %>%
    as.data.frame()
  colnames(species.lvl) <- species.lvl["species",]
  species.lvl <- species.lvl[rownames(species.lvl)!="species",]
  
  # get host percent
  genus.lvl$`Percent host` <- get_percent_host(rownames(genus.lvl))
  species.lvl$`Percent host` <- get_percent_host(rownames(species.lvl))
  
  # plot and save genus level
  genus.plot <- plot_data(genus.lvl)
  
  save.path <- file.path(opts$output_dir, 
                         paste0(opts$taxonomy_method, "_",
                                "genus.pdf"))
  ggsave(save.path, genus.plot)
  
  
  # plot and save species level
  species.plot <- plot_data(species.lvl)
  save.path <- file.path(opts$output_dir, 
                         paste0(opts$taxonomy_method, "_",
                                "species.pdf"))
  ggsave(save.path, species.plot)
  
}

main()

