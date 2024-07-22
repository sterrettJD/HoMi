library(tidyverse)
library(ggplot2)
library(ggbeeswarm)

df <- read.csv("reads_breakdown.csv")
df <- dplyr::rename(df, 
              sample_name=X)

df$true_perc_host <- df$sample_name %>% 
  gsub(pattern="_perc_host_\\d", replacement="") %>%
  as.numeric()

df$Percent.host <- df$Percent.host*100

mod <- lm(Percent.host ~ true_perc_host, data=df)
print(summary(mod))

ggplot(df, mapping=aes(x=true_perc_host, y=Percent.host)) +
  geom_jitter(width=1, size=3, alpha=0.8) +
  geom_smooth(method="lm") +
  theme_bw(base_size=22) +
  labs(x="True percent host reads (jittered)", y="Actual percent host reads")

ggsave("percent_host_reads_plot.pdf")
