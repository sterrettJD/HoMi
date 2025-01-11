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
    make_option(c("-m", "--metadata_file"),
                help="The metadata file for HoMi",
                default=NULL),
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
}

main()

