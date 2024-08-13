library(tidyverse)
library(ggplot2)
library(ggbeeswarm)


get_args <- function(){
  option_list <- list( 
      make_option(c("-i", "--input_file"),
                  help="The reads breakdown CSV output by HoMi"),
      make_option(c("-o", "--output_plot"), 
                  help="The path to create the output plot")

  )
  opt <- parse_args(OptionParser(option_list=option_list))
  return (opt)
}


clean_df <- function(df){
  df <- dplyr::rename(df, 
                      sample_name=X)

  df$true_perc_host <- df$sample_name %>% 
    gsub(pattern="_perc_host_\\d", replacement="") %>%
    as.numeric()

  df$Percent.host <- df$Percent.host*100

  return (df)
}


main <- function(){
  args <- get_args()
  df <- read.csv(args$input_file)
  df <- clean_df
  
  mod <- lm(Percent.host ~ true_perc_host, data=df)
  print(summary(mod))

  ggplot(df, mapping=aes(x=true_perc_host, y=Percent.host)) +
    geom_jitter(width=1, size=3, alpha=0.8) +
    geom_smooth(method="lm") +
    theme_bw(base_size=22) +
    labs(x="True percent host reads (jittered)", y="Actual percent host reads")

  ggsave(args$output_plot)
}


main()