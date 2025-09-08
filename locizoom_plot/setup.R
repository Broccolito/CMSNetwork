suppressPackageStartupMessages({
  if(!require("dplyr")){
    install.packages("dplyr")
    library("dplyr")
  }
  
  if(!require("ggplot2")){
    install.packages("ggplot2")
    library("ggplot2")
  }
  
  if(!require("ggpubr")){
    install.packages("ggpubr")
    library("ggpubr")
  }
  
  if(!require("latex2exp")){
    install.packages("latex2exp")
    library("latex2exp")
  }
})

if(!exists("cms_hg19")){
  cat("Loading hg19 database...\n")
  load("cms_hg19.rda")
}
