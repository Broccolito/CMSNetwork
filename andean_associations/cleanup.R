library(dplyr)
library(vcfR)
library(adegenet)
library(data.table)
library(readxl)
library(writexl)
library(purrr)

vcf = read.vcfR("Andean_CDP_QC.vcf.gz")
gen = vcfR2genlight(vcf)
gen = as.data.frame(gen)
metadata = fread("metadata.csv")
overlap = read_excel(path = "tibetan_andean_overlap_significant.xlsx")
wgs_db = fread("wgs_db.csv")

pheno = read_excel(path = "wgs_pheno.xlsx")
hvr_data = fread("hvr_data.csv")
sleep_data = fread("sleep_data.csv")

variant_ids = vector()
variant_lookup = vector()
for(i in 1:nrow(overlap)){
  variant_chr = overlap$chr_andean[i]
  variant_pos = overlap$pos_hg37_andean[i]
  vcf_index = which(vcf@fix[,1] == variant_chr & vcf@fix[,2] == variant_pos)
  variant_id = vcf@fix[vcf_index,3]
  
  variant_lookup = rbind.data.frame(variant_lookup,
                                    data.frame(
                                      chr = variant_chr, 
                                      pos = variant_pos,
                                      variant_id = variant_id,
                                      gene = overlap$gene_symbol[i]
                                    )
  )
  
  variant_ids = c(variant_ids, variant_id)
}
variant_ids = unique(variant_ids)

write.csv(variant_lookup, file = "variant_lookup.csv", 
          row.names = FALSE, quote = FALSE)

genotypes = gen[, names(gen) %in% variant_ids]
wgs_id = rownames(genotypes)
wgs_id = strsplit(wgs_id, split = "_") |>
  lapply(function(x){x[1]}) |>
  unlist()
genotypes$wgs_id = wgs_id
genotypes = inner_join(genotypes, metadata, by = "wgs_id")

genotypes = inner_join(genotypes, select(wgs_db, iid, phv00001), 
                       by = join_by("initials" == "phv00001"))

geno_pheno_master = inner_join(genotypes, pheno, by = "iid")
geno_pheno_hvr = inner_join(genotypes, hvr_data,  by = "iid")
geno_pheno_sleep = inner_join(genotypes, sleep_data, by = "iid")

write.csv(geno_pheno_master, file = "geno_pheno_master_uncleaned.csv")
write.csv(geno_pheno_hvr, file = "geno_pheno_hvr_uncleaned.csv")
write.csv(geno_pheno_sleep, file = "geno_pheno_sleep_uncleaned.csv")


