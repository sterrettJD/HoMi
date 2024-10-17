pseq_from_bracken <- function(k, as.relative=TRUE){
    k.tax.table <- k[c("domain", "phylum", "class", 
                       "order", "family", "genus", "species")]
    
    if(as.relative==T){
        rel.cols <- colnames(k)[grepl(x=colnames(k), pattern="_frac")]
        k.feature.table <- k[,rel.cols]
        colnames(k.feature.table) <- gsub(
            colnames(k.feature.table), 
            pattern=".bracken_frac", replacement="")
    } 
    else if (as.relative==F) {
        num.cols <- colnames(k)[grepl(x=colnames(k), pattern="_frac")]
        k.feature.table <- k[,num.cols]
        colnames(k.feature.table) <- gsub(
            colnames(k.feature.table), 
            pattern=".bracken_num", replacement="")
    }
    k.pseq <- phyloseq::phyloseq(phyloseq::otu_table(k.feature.table, 
                                                     taxa_are_rows=T),
                       phyloseq::tax_table(as.matrix(k.tax.table)))

    return(k.pseq)
}