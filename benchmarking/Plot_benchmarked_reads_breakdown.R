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
                  help="The reads breakdown CSV output by HoMi"),
      make_option(c("-m", "--metadata_file"),
                  help="The metadata file for HoMi",
                  default=NULL),
      make_option(c("-c", "--column_to_use"),
                  help="The column to use in the plot",
                  default="X"),
      make_option(c("-n", "--name_for_plot"),
                  help="What to title the specified column in the output plot",
                  default="True percent host reads (jittered)"),
      make_option(c("-j", "--jitter_width"),
                  help="The width to jitter the x axis. Should be an integer.",
                  default=1, type="integer"),
      make_option(c("--no_dotted_line"),
                  help="Pass this to avoid showing the dotted line on the benchmark plot.",
                  default=FALSE, action="store_true"),
      make_option(c("-o", "--output_plot"), 
                  help="The path to create the output plot")

  )
  opt <- parse_args(OptionParser(option_list=option_list))
  return (opt)
}


clean_df <- function(df, col_to_use=X){
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
  df <- read.csv(args$input_file)
  if (!is.null(args$metadata_file)){
    metadata <- read.csv(args$metadata_file)
    df <- merge(df, metadata, by.x="X", by.y="Sample")
  }
  
  df <- clean_df(df, col_to_use=args$column_to_use)

  mod <- lm(Percent.host ~ true_perc_host, data=df)
  print(summary(mod))

  p <- ggplot(df, mapping=aes(x=true_perc_host, y=Percent.host)) +
    geom_jitter(width=args$jitter_width, size=3, alpha=0.8) +
    geom_smooth(method="lm") +
    theme_bw(base_size=22) +
    xlim(-5, 90) +
    ylim(-5, 90) +
    labs(x=args$name_for_plot, y="Recovered percent host reads")

  if (!args$no_dotted_line) {
    p <- p + 
      geom_abline(slope=1, intercept=0, 
      linetype="dotted", color="black", size=1)
  }

  ggsave(args$output_plot, p)
}


main()
