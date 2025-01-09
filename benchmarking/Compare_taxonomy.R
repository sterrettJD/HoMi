if (!require("tidyverse")){
  install.packages("tidyverse", repos="http://cran.us.r-project.org")
  library("tidyverse")
}

if (!require("optparse")){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library("optparse")
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
    select(all_of(tax.columns)) %>%
    as.matrix()
  rownames(taxonomy.table) <- bugs
  
  
  return(list(abundance=feature.table, tax=taxonomy.table))
}


main <- function(){
  opts <- get_args()
  data <- read_taxonomy_data(opts$input_file, type=opts$taxonomy_method)  
  str(data$abundance)
  print(rownames(data$abundance)[1:10])
  str(data$tax)
}

main()

