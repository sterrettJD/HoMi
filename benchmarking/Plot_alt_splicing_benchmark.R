if (!require("tidyverse")){
  install.packages("tidyverse", repos="http://cran.us.r-project.org")
  library("tidyverse")
}

if (!require("ggplot2")){
  install.packages("ggplot2", repos="http://cran.us.r-project.org")
  library("ggplot2")
}

if (!require("optparse")){
  install.packages("optparse", repos="http://cran.us.r-project.org")
  library("optparse")
}


get_args <- function(){
  option_list <- list( 
      make_option(c("-i", "--input_file"),
                  help="The aggregated hostile results file"),
      make_option(c("-c", "--column_of_sample_names"),
                  help="The column to use with sample names for percent host parsing",
                  default="fastq1_in_name"),
      make_option(c("-n", "--name_for_plot"),
                  help="What to title the specified column in the output plot",
                  default="True percent host reads (jittered)"),
      make_option(c("-j", "--jitter_width"),
                  help="The width to jitter the x axis. Should be an integer.",
                  default=1, type="integer"),
      make_option(c("--no_dotted_line"),
                  help="Pass this to avoid showing the dotted line on the benchmark plot.",
                  default=FALSE, action="store_true"),
      make_option(c("-l", "--output_lms"), 
                  help="The path to create the output lm results csv"),
      make_option(c("-o", "--output_plot"), 
                  help="The path to create the output plot")

  )
  opt <- parse_args(OptionParser(option_list=option_list))
  return (opt)
}

get_simulation_method <- function(filepath){
  patterns <- c("synthetic_communities", "synthetic_transcriptomes", "semi")
  for(pattern in patterns){
    if(grepl(pattern, filepath)){
      return(pattern)
    }   
  }
  stop("None of the accepted patterns were detected in ", filepath)
}


clean_df <- function(df, col_to_use=X){
  df <- dplyr::rename(df, 
                      sample_name=col_to_use)
  colnames(df)[1] <- "filepath"
  df$simulation_method <- sapply(df$filepath, get_simulation_method)
  
  df[df$simulation_method=="semi","simulation_method"] <- "semisynthetic_transcriptomes"
  
  
  df$true_perc_host <- df$sample_name %>% 
    str_split_i(pattern="_",i=1) %>%
    as.numeric()

  df$Percent_host <- df$reads_removed_proportion*100
  
  return (df)
}


main <- function(){
  args <- get_args()
  df <- read.csv(args$input_file)
  df <- clean_df(df, col_to_use=args$column_of_sample_names)

  run_regression <- function(df){broom::tidy(summary(lm(Percent_host ~ true_perc_host, data=df)))}
  reg_summaries <- plyr::ddply(df, ~aligner+simulation_method, run_regression)
  print(reg_summaries)
  write.csv(x=reg_summaries, file=args$output_lms)

  p <- ggplot(df, mapping=aes(x=true_perc_host, y=Percent_host, color=aligner)) +
    geom_jitter(width=args$jitter_width, size=3, alpha=0.8) +
    geom_smooth(method="lm") +
    theme_bw(base_size=22) +
    xlim(-5, 90) +
    ylim(-5, 90) +
    facet_grid(cols=vars(simulation_method))
    labs(x=args$name_for_plot, y="Recovered percent host reads")

  if (!args$no_dotted_line) {
    p <- p + 
      geom_abline(slope=1, intercept=0, 
      linetype="dotted", color="black", size=1)
  }

  ggsave(args$output_plot, p)
}


main()
