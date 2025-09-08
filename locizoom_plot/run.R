get_directory = function(){
  args = commandArgs(trailingOnly = FALSE)
  file = "--file="
  rstudio = "RStudio"
  match = grep(rstudio, args)
  if(length(match) > 0){
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }else{
    match = grep(file, args)
    if (length(match) > 0) {
      return(dirname(normalizePath(sub(file, "", args[match]))))
    }else{
      return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
    }
  }
}

setwd(get_directory())
rm(list = ls())
graphics.off()
gc()
source("setup.R")
source("locizoom.R")

locizoom(gene_name = "",                    # Gene name
         rs_num = "",                       # RS number and protein variant (if applicable)
         snp_chr = 2,                       # Chromosome
         snp_loc = 46588031,                # Location (Build hg19)
         threshold = 3e5,                   # Window size
         file_name = "EPAS1 Locizoom.png",  # Saved file name
         plot_ratio = 0.3,                  # Locizoom plot / gene map ratio
         plot_height = 5,                   # Inch
         plot_width = 7,                    # Inch
         text_size = 15,                    # Font size in MS word units
         text_shift = 4,                    # How much the gene annotation is shifted up
         text_horizontal_shift = 0,         # How much the gene annotation is shifted horizontally
         include_genetrack = FALSE,         # Whether to include the gene tracks
         include_tibetan = TRUE,            # Whether to include the Tibetans
         andean_color = "#D39200",          # Color of the Andean markers
         tibetan_color = "#7CAE00",         # Color of the Tibetan markers
         save_plot = TRUE                   # Save plot
) |>
  write.csv(file = "saved_markers/markers.csv", 
            quote = FALSE, row.names = FALSE)

