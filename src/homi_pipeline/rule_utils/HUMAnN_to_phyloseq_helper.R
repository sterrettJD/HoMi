if(!require("stringi")){
    install.packages("stringi", repos="http://cran.us.r-project.org")
}
if(!require("r2r")){
    install.packages("r2r", repos="http://cran.us.r-project.org")
}


# These should already be installed 
library(data.table)
library(tidyverse)
library(stringi)
library(r2r)



filter_only_genes <- function(df, with_taxa=F){
    feature.names <- rownames(df)
    # if there's a | in the rowname, it contains a corresponding species
    contains.taxa <- grepl(x=feature.names, pattern="\\|")
    
    if (with_taxa){
        to.return <- gn[contains.taxa,]
        rownames(to.return) <- feature.names[contains.taxa]
    } else {
        to.return <- gn[!contains.taxa,]
        rownames(to.return) <- feature.names[!contains.taxa]
    }
    
    return(to.return)
}


get_ec_class_db <- function(){
    if(!file.exists("enzclass.txt")){
        download.file("https://ftp.expasy.org/databases/enzyme/enzclass.txt",
                      destfile="enzclass.txt")
        
    } 
    
    enzclass <- read_file("enzclass.txt")
    # split by lines and remove header/footer
    splitted <- str_split(enzclass, pattern="\n")[[1]][12:418]
    # split EC from the corresponding name
    parsed <- sapply(splitted, 
           FUN=function(x) str_split(x, pattern="  ") %>% 
               unlist() %>%
               stri_remove_empty()
           )
    df <- matrix(parsed, nrow=2) %>% t() %>% as.data.frame()
    colnames(df) <- c("EC", "name")
    # trim ECs for hierarchical mapping
    df$EC.trimmed <- df$EC %>% 
        sapply(FUN=function(x) str_split(x, pattern="-")[[1]][1] %>%
                   gsub(pattern="\\ ", replacement=""))
    
    # create a hashmap to map keys and values
    mapper <- hashmap()
    mapper[df$EC.trimmed] <- df$name
    return(mapper)
}

extract_ecs <- function(genes){
    # in the case of multiple ECs, this only extracts the first. 
    # Should be good enough for general plotting though.
    # Regex grabs the pattern num.num.num.num, like 1.1.1.1
    return(str_extract(genes, 
                        pattern="([1-9|-]+\\.){3}[1-9|-]+"))
}

get_ecs_hierarchy <- function(ecs){
    cbind(str_extract(ecs, 
                pattern="([1-9|-]+\\.)"),
          str_extract(ecs, 
                      pattern="([1-9|-]+\\.){2}"),
        
          str_extract(ecs, 
                      pattern="([1-9|-]+\\.){3}")
          )
}

make_pseq_tax_table_from_genefams <- function(genes){
    mapper <- get_ec_class_db()
    
    ecs <- extract_ecs(genes)
    ecs.hier <- get_ecs_hierarchy(ecs)
    
    ecs.hier.names <- data.frame(matrix(nrow=length(genes)))
    ecs.hier.names$Level.1 <- mapper[ecs.hier[,1]]
    ecs.hier.names$Level.2 <- mapper[ecs.hier[,2]]
    ecs.hier.names$Level.3 <- mapper[ecs.hier[,3]]
    ecs.hier.names$matrix.nrow...length.genes.. <- NULL
    
    cleaned.names <- sapply(genes, 
                           FUN=function(x) str_split(x, pattern="\\(expasy\\) ")[[1]][2])
    
    ecs.hier.names$Level.4 <- cleaned.names
    rownames(ecs.hier.names) <- genes
    
    return(ecs.hier.names)
}
