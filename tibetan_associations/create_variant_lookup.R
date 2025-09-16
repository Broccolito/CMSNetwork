library(dplyr)
library(data.table)
library(readxl)
library(writexl)
library(vcfR)
library(adegenet)

overlap_variants = read_excel(path = "tibetan_andean_overlap_significant.xlsx")

file_list = list.files(path = "vcf", pattern = "vcf.gz", full.names = TRUE)
file_list = file_list[!grepl(pattern = ".tbi", file_list)]

vcf = list()
for(f in file_list){
  chr_name = strsplit(f, split = "[_]") |>
    lapply(function(x){x[2]}) |>
    unlist()
  cat(paste0("Reading genetics data from ", chr_name, "...\n"))
  vcf[[chr_name]] = read.vcfR(f)
}

gen = list()
for(f in file_list){
  chr_name = strsplit(f, split = "[_]") |>
    lapply(function(x){x[2]}) |>
    unlist()
  cat(paste0("Reading genetics data from ", chr_name, " into genlight format...\n"))
  gen[[chr_name]] = vcfR2genlight(read.vcfR(f)) |>
    as.data.frame()
}

vcf = vcf[paste0("chr", 1:22)]
gen = gen[paste0("chr", 1:22)]


variant_lookup = vector()

for(i in 1:nrow(overlap_variants)){
  chr = overlap_variants$chr_tibetan[i]
  chr = paste0("chr", chr)
  pos = overlap_variants$pos_hg37_tibetan[i]
  gene = overlap_variants$gene_symbol[i]
  fix = vcf[[chr]]@fix
  
  if(any(fix[,2] == as.character(pos))){
    variant_index = which(fix[,2] == as.character(pos))
    genotype = gen[[chr]][,variant_index]
    selected_variant_r2 = 1
    selected_variant_gene = gene
    selected_variant_pos = pos
    
  }else{
    
    haplotype_annotation = file.path("variant_annotation", paste0(tolower(gene), "_variant_annotation.txt")) |>
      fread() |>
      as.data.frame()
    
    haplotype_annotation = haplotype_annotation |>
      arrange(desc(R2)) |>
      distinct(Position, .keep_all = TRUE) |>
      # filter(R2 >= 0.6) |>
      mutate(pos = gsub(pattern = paste0(chr, ":"), replacement = "", Position))
    
    haplotype_overlap = data.frame(pos = fix[,2]) |>
      inner_join(haplotype_annotation, by = "pos") |>
      select(pos, R2, `Gene Symbol`)
    
    if(nrow(haplotype_overlap) > 0){
      
      if(haplotype_overlap$R2[1] >= 0.7){
        selected_variant_pos = haplotype_overlap$pos[1]
        selected_variant_r2 = haplotype_overlap$R2[1]
        selected_variant_gene = haplotype_overlap[["Gene Symbol"]][1]
      }else{
        if(gene %in% haplotype_overlap[["Gene Symbol"]]){
          haplotype_overlap = haplotype_overlap |>
            filter(`Gene Symbol` == gene)
        }
        selected_variant_pos = haplotype_overlap$pos[1]
        selected_variant_r2 = haplotype_overlap$R2[1]
        selected_variant_gene = haplotype_overlap[["Gene Symbol"]][1]
      }
      
    }else{
      selected_variant_r2 = NA
      selected_variant_gene = NA
      selected_variant_pos = NA
    }
    
  }
  
  if(nrow(haplotype_overlap) > 0){
    variant_index = which(fix[,2] %in% selected_variant_pos)[1]
    variant_lookup = rbind.data.frame(variant_lookup,
                                      data.frame(
                                        chr = chr,
                                        pos = pos,
                                        gene = gene,
                                        selected_variant_pos = selected_variant_pos,
                                        selected_variant_r2 = selected_variant_r2,
                                        selected_variant_gene = selected_variant_gene
                                      )
    )
  }
  
  rm(selected_variant_r2)
  rm(selected_variant_gene)
  rm(selected_variant_pos)
  
}

write_xlsx(variant_lookup, path = "tibetan_variant_lookup.xlsx")




