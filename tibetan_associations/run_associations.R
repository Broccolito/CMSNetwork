library(dplyr)
library(data.table)
library(readxl)
library(writexl)
library(vcfR)
library(adegenet)

overlap_variants = read_excel(path = "tibetan_andean_overlap_significant.xlsx")
variant_lookup = read_excel(path = "tibetan_variant_lookup.xlsx")
phenotypes = read_excel("nepali_physiology.xlsx", sheet = "Sorted by Menstruation")

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

genotypes = list()
for(i in 1:nrow(variant_lookup)){
  chr = variant_lookup$chr[i]
  pos = variant_lookup$selected_variant_pos[i]
  gene = variant_lookup$gene[i]
  genotype = gen[[chr]][, which(vcf[[chr]]@fix[,2] == pos)]
  genotypes[[gene]] = genotype
}

genotypes = do.call(cbind.data.frame, genotypes)
subject_id = colnames(vcf[["chr1"]]@gt)[-1]
subject_id = gsub(pattern = "TIBETN_ADR", replacement = "", subject_id) |>
  as.numeric()
genotypes[["subject_id"]] = subject_id
genotypes = genotypes |>
  select(subject_id, everything())

names(phenotypes)[3] = "subject_id" 
master = inner_join(phenotypes, genotypes, by = "subject_id")

phenotype_names = select(master, `Inspired O2`:`AVG HVR BTPS`, `AVG HHRR...58`, height:BMI) |>
  names()
genotype_names = select(master, CTSB:TBC1D7) |>
  names()

master_regression_stats = vector()
for(g in genotype_names){
  cat(paste0("\n\n\nRegressing phenotypes on the variant in the ", g, " gene...\n\n\n"))
  for(p in phenotype_names){
    cat(paste0("Testing associations with ", p, "...\n"))
    
    try({
      
      regression_data = data.frame(
        geno = master[[g]],
        pheno = master[[p]],
        age = master$age,
        menstration = master[["Still Menstruating?"]]
      )
      
      data_premenopause = regression_data |>
        filter(menstration == "YES")
      data_postmenopause = regression_data |>
        filter(menstration == "NO")
      
      n_total = nrow(regression_data)
      suppressWarnings({
        regression_data = regression_data |>
          mutate(phenotype = as.numeric(pheno)) |>
          mutate(age = as.numeric(age)) |>
          na.omit()
      })
      n_avail = nrow(regression_data)
      n0 = sum(regression_data$geno == 0)
      n1 = sum(regression_data$geno == 1)
      n2 = sum(regression_data$geno == 2)
      hetero = length(unique(regression_data$geno)) > 1
      
      l = lm(pheno ~ geno + age + menstration, regression_data)
      l = summary(l)
      age_beta = l$coefficients[3,1]
      age_p = l$coefficients[3,4]
      variant_beta = l$coefficients[2,1]
      variant_p = l$coefficients[2,4]
      menstration_beta = l$coefficients[4,1]
      menstration_p = l$coefficients[4,4]
      
      l_pre = lm(pheno ~ geno + age, data_premenopause)
      l_pre = summary(l_pre)
      pre_age_beta = l_pre$coefficients[3,1]
      pre_age_p = l_pre$coefficients[3,4]
      pre_variant_beta = l_pre$coefficients[2,1]
      pre_variant_p = l_pre$coefficients[2,4]
      
      l_post = lm(pheno ~ geno + age, data_postmenopause)
      l_post = summary(l_post)
      post_age_beta = l_post$coefficients[3,1]
      post_age_p = l_post$coefficients[3,4]
      post_variant_beta = l_post$coefficients[2,1]
      post_variant_p = l_post$coefficients[2,4]
      
      master_regression_stat = data.frame(
        gene = g,
        phenotype = p,
        n_total = n_total,
        n_avail = n_avail,
        n0, n1, n2,
        variant_beta,
        variant_p, 
        pre_variant_beta,
        pre_variant_p,
        post_variant_beta,
        post_variant_p,
        age_beta,
        age_p,
        menstration_beta,
        menstration_p,
        pre_age_beta,
        pre_age_p,
        post_age_beta,
        post_age_p
      )
      
      master_regression_stats = rbind.data.frame(
        master_regression_stats,
        master_regression_stat
      )
      
    }, silent = TRUE)
    
  }
}

master_regression_stats = master_regression_stats |>
  arrange(variant_p)

sig_regression_stats = master_regression_stats |>
  filter(pre_variant_p <= 0.05 | post_variant_p <= 0.05 | variant_p <= 0.05) |>
  arrange(gene, variant_p) |>
  filter(!(phenotype %in% c("Body Temp (°K)", "Amb Temp (°K)", "Humidity (fraction)",
                            "PH2O at amb T (mmHg)", "PH2O at Body Temp (mmHg)", 
                            "Barometric Pressure (mmHg)", "3L Cal Syringe value")))

write_xlsx(sig_regression_stats, path = "sig_association_stats.xlsx")
write_xlsx(master_regression_stats, path = "all_association_stats.xlsx")

