if (!require("tidyverse")) {
    install.packages("tidyverse", repos="http://cran.us.r-project.org")
}

if (!require("data.table")) {
    install.packages("data.table", repos="http://cran.us.r-project.org")
}

if (!require("ggplot2")) {
    install.packages("ggplot2", repos="http://cran.us.r-project.org")
}

if (!require("patchwork")) {
    install.packages("patchwork", repos="http://cran.us.r-project.org")
}

if (!require("cowplot")) {
    install.packages("cowplot", repos="http://cran.us.r-project.org")
}

if (!require("remotes")) {
    install.packages("remotes", repos="http://cran.us.r-project.org")
}

if (!require("phyloseq")) {
    BiocManager::install("phyloseq")
}

if (!require("microshades")) {
    remotes::install_github("KarstensLab/microshades")
}

if (!require("speedyseq")) {
    remotes::install_github("mikemc/speedyseq")
}

if (!require("Nonpareil")){
  install.packages('Nonpareil', repos="http://cran.us.r-project.org")
}

file.create("R_packages_installed")
