if (!require("tidyverse")){
  install.packages("tidyverse", repos="http://cran.us.r-project.org")
  library("tidyverse")
}

if (!require("ggplot2")){
  install.packages("ggplot2", repos="http://cran.us.r-project.org")
  library("ggplot2")
}

if (!require("ggbeeswarm")){
  install.packages("ggbeeswarm", repos="http://cran.us.r-project.org")
  library("ggbeeswarm")
}

if (!require("optparse")){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library("optparse")
}


get_args <- function(){
  option_list <- list( 
    make_option(c("-i", "--input_file"),
                help=paste0("A CSV file with the following columns: ",
                            "input_file, hostile_index, dataset, ",
                            "metadata_file, column_to_use, name_for_plot, ",
                            "jitter_width, no_dotted_line. ",
                            "Column descriptions can be found in", 
                            "`Plot_benchmarked_reads_breakdown.R` arguments.")),
    make_option(c("-o", "--output_plot"), 
                help="The path to create the output plot")
    
  )
  opt <- parse_args(OptionParser(option_list=option_list))
  return (opt)
}


clean_df <- function(df, col_to_use="X"){
  df <- dplyr::rename(df, 
                      sample_name=col_to_use)
  
  df$true_perc_host <- df$sample_name %>% 
    gsub(pattern="_perc_host_\\d", replacement="") %>%
    as.numeric()
  
  df$Percent.host <- df$Percent.host*100
  
  return (df)
}


main <- function(){
  args <- get_args()
  meta.df <- read.csv(args$input_file)
  
  big.data <- data.frame(matrix(ncol=11, nrow=0))
  for (project in meta.df["input_file"]){
    # read and process reads breakdown file
    notclean.df <- read.csv(project)
    clean.df <- clean_df(notclean.df)
    
    # add project name to df
    clean.df$input_file <- project
    clean.df$dataset <- meta.df[meta.df["input_file"]==project, 
                                "dataset"]
    clean.df$hostile_index <- meta.df[meta.df["input_file"]==project, 
                                      "hostile_index"]
    
    dataset.metadata.file <- meta.df[meta.df["input_file"]==project, "metadata_file"]
    column.to.use <- meta.df[meta.df["input_file"]==project, "column_to_use"]
    
    if (!is.na(column.to.use)){
      dataset.metadata <- read.csv(dataset.metadata.file)
      # replace percent host with better percent host
      # by merging this and then renaming
      clean.df <- merge(clean.df, dataset.metadata[c("Sample", column.to.use)],
                      by.x="sample_name", by.y="Sample")
      clean.df$plotting_true_percent_host <- clean.df[column.to.use]
    } else {
      clean.df$plotting_true_percent_host <- clean.df$true_percent_host
    }
    
    # rbind with the full dataframe
    big.data <- rbind(big.data, clean.df)
  }
  
  
  # Get list of facet values to draw theoretical optimal line on
  facet_vals_with_line <- unique(big.data$dataset)
  facet_vals_with_line <- facet_vals_with_line[facet_vals_with_line!="Pereira-Marques"]  # replace with actual name
  
  # Create a small data frame to use with geom_abline
  abline_df <- data.frame(dataset=facet_vals_with_line)
  
  p <- ggplot(big.data, 
              mapping=aes(x=plotting_true_percent_host, 
                          y=Percent.host, 
                          color=hostile_index)) +
    geom_jitter(width=args$jitter_width, size=3, alpha=0.8) +
    geom_smooth(method="lm") +
    theme_bw(base_size=22) +
    xlim(-5, 90) +
    ylim(-5, 90) +
    labs(x="True percent host reads", y="Recovered percent host reads") +
    facet_grid(cols=dataset, 
               labeller=label_wrap_gen(10)) +
    geom_abline(data=abline_df, 
                aes(slope=1, intercept=0), 
                linetype="dashed", 
                color="black", 
                inherit.aes=FALSE)
   
}


main()
