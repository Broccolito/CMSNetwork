library(dplyr)
library(vcfR)
library(data.table)
library(purrr)
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)

variant_lookup = fread("variant_lookup.csv")
pheno_hvr = fread("geno_pheno_hvr_cleaned.csv")
pheno_sleep = fread("geno_pheno_sleep_cleaned.csv")
pheno_master = fread("geno_pheno_master_cleaned.csv")
associations = read_excel(path = "Andean Associations.xlsx")
vcf = read.vcfR("Andean_CDP_QC.vcf.gz") 

for(i in 1:nrow(associations)){
  
  linkage = associations[i,]
  
  variant_name = linkage$variant
  phenotype_name = linkage$phenotype
  pheno_data_name = linkage$table
  
  if(pheno_data_name == "master"){
    pheno_data = pheno_master
  }else{
    if(pheno_data_name == "hvr"){
      pheno_data = pheno_hvr
    }else{
      pheno_data = pheno_sleep
    }
  }
  
  data = data.frame(
    phenotype = pheno_data[[phenotype_name]],
    variant = pheno_data[[variant_name]],
    sex = pheno_data$sex,
    age = pheno_data$age
  ) |>
    na.omit()
  
  variant_info = vcf@fix[which(vcf@fix[,3] == variant_name),]
  ref = variant_info[4]
  alt = variant_info[5]
  
  plt = ggplot(data, aes(x = as.factor(variant), y = phenotype)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(fill = sex), color = "black", shape = 21, size = 3, 
               position = position_dodge(0.4), alpha = 0.7) + 
    scale_fill_manual(values = c("#B22222", "#000080")) + 
    ggtitle(label = NULL, subtitle = paste0("p = ", round(linkage$variant_p, 3))) + 
    xlab(paste0(linkage$gene, " ", variant_name, " (", ref, ">", alt, ")")) + 
    ylab(phenotype_name) + 
    theme_classic2() + 
    theme(text = element_text(size = 15), legend.position = "none"); plt
  
  plt_name = paste0(phenotype_name, " vs. ", linkage$gene, ".png")
  ggsave(filename = file.path("figures", plt_name), plot = plt, device = "png", dpi = 1200,
         width = 3.5, height = 4.5)
  
}

