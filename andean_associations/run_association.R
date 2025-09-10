library(dplyr)
library(data.table)
library(purrr)

variant_lookup = fread("variant_lookup.csv")
pheno_hvr = fread("geno_pheno_hvr_cleaned.csv")
pheno_sleep = fread("geno_pheno_sleep_cleaned.csv")
pheno_master = fread("geno_pheno_master_cleaned.csv")

hvr_phenotypes = select(pheno_hvr, height:abpm_conventional_dbp) |>
  names()

sleep_phenotypes = select(pheno_sleep, height:hrr_co2) |>
  names()

master_phenotypes = select(pheno_master, height_meters:free_testosterone) |>
  names()

run_regression_by_variant = function(variant_id){
  
  master_regression_stats = vector()
  for(pheno in master_phenotypes){
    regression_data = data.frame(
      variant = pheno_master[[variant_id]],
      sex = pheno_master[["sex"]],
      age = pheno_master[["age"]],
      phenotype = pheno_master[[pheno]]
    ) 
    
    n_total = nrow(regression_data)
    suppressWarnings({
      regression_data = regression_data |>
        mutate(phenotype = as.numeric(phenotype)) |>
        na.omit()
    })
    n_avail = nrow(regression_data)
    n0 = sum(regression_data$variant == 0)
    n1 = sum(regression_data$variant == 1)
    n2 = sum(regression_data$variant == 2)
    hetero = length(unique(regression_data$variant)) > 1
    
    if(n_avail > 3 & hetero){
      if(all(regression_data$sex == "M") | all(regression_data$sex == "F")){
        l = lm(phenotype ~ variant + age, regression_data)
        l = summary(l)
        sex_beta = NA
        sex_p = NA
        age_beta = l$coefficients[3,1]
        age_p = l$coefficients[3,4]
      }else{
        l = lm(phenotype ~ variant + sex + age, regression_data)
        l = summary(l)
        sex_beta = l$coefficients[3,1]
        sex_p = l$coefficients[3,4]
        age_beta = l$coefficients[4,1]
        age_p = l$coefficients[4,4]
      }
      
      nm = sum(regression_data$sex == "M")
      nf = sum(regression_data$sex == "F")
      
      variant_beta = l$coefficients[2,1]
      variant_p = l$coefficients[2,4]
      
      master_regression_stat = data.frame(
        variant = variant_id,
        phenotype = pheno,
        n_total = n_total,
        n_avail = n_avail,
        nm = nm, 
        nf = nf,
        n0 = n0,
        n1 = n1,
        n2 = n2,
        variant_beta = variant_beta,
        variant_p = variant_p,
        sex_beta = sex_beta,
        sex_p = sex_p,
        age_beta = age_beta,
        age_p = age_p
      )
      
      master_regression_stats = rbind.data.frame(
        master_regression_stats,
        master_regression_stat
      )
    }else{
      cat(paste0("Insufficiant sample size for ", pheno, " to conduct a linear regression...\n"))
    }
    
  }
  
  sleep_regression_stats = vector()
  for(pheno in sleep_phenotypes){
    regression_data = data.frame(
      variant = pheno_sleep[[variant_id]],
      sex = pheno_sleep[["sex"]],
      age = pheno_sleep[["age"]],
      phenotype = pheno_sleep[[pheno]]
    ) 
    
    n_total = nrow(regression_data)
    suppressWarnings({
      regression_data = regression_data |>
        mutate(phenotype = as.numeric(phenotype)) |>
        na.omit()
    })
    n_avail = nrow(regression_data)
    
    n0 = sum(regression_data$variant == 0)
    n1 = sum(regression_data$variant == 1)
    n2 = sum(regression_data$variant == 2)
    hetero = length(unique(regression_data$variant)) > 1
    
    if(n_avail > 3 & hetero){
      if(all(regression_data$sex == "M") | all(regression_data$sex == "F")){
        l = lm(phenotype ~ variant + age, regression_data)
        l = summary(l)
        sex_beta = NA
        sex_p = NA
        age_beta = l$coefficients[3,1]
        age_p = l$coefficients[3,4]
      }else{
        l = lm(phenotype ~ variant + sex + age, regression_data)
        l = summary(l)
        sex_beta = l$coefficients[3,1]
        sex_p = l$coefficients[3,4]
        age_beta = l$coefficients[4,1]
        age_p = l$coefficients[4,4]
      }
      
      nm = sum(regression_data$sex == "M")
      nf = sum(regression_data$sex == "F")
      
      variant_beta = l$coefficients[2,1]
      variant_p = l$coefficients[2,4]
      
      sleep_regression_stat = data.frame(
        variant = variant_id,
        phenotype = pheno,
        n_total = n_total,
        n_avail = n_avail,
        nm = nm, 
        nf = nf,
        n0 = n0,
        n1 = n1,
        n2 = n2,
        variant_beta = variant_beta,
        variant_p = variant_p,
        sex_beta = sex_beta,
        sex_p = sex_p,
        age_beta = age_beta,
        age_p = age_p
      )
      
      sleep_regression_stats = rbind.data.frame(
        sleep_regression_stats,
        sleep_regression_stat
      )
    }else{
      cat(paste0("Insufficiant sample size for ", pheno, " to conduct a linear regression...\n"))
    }
    
  }
  
  hvr_regression_stats = vector()
  for(pheno in hvr_phenotypes){
    regression_data = data.frame(
      variant = pheno_hvr[[variant_id]],
      sex = pheno_hvr[["sex"]],
      age = pheno_hvr[["age"]],
      phenotype = pheno_hvr[[pheno]]
    ) 
    
    n_total = nrow(regression_data)
    suppressWarnings({
      regression_data = regression_data |>
        mutate(phenotype = as.numeric(phenotype)) |>
        na.omit()
    })
    n_avail = nrow(regression_data)
    n0 = sum(regression_data$variant == 0)
    n1 = sum(regression_data$variant == 1)
    n2 = sum(regression_data$variant == 2)
    hetero = length(unique(regression_data$variant)) > 1
    
    if(n_avail > 3 & hetero){
      if(all(regression_data$sex == "M") | all(regression_data$sex == "F")){
        l = lm(phenotype ~ variant + age, regression_data)
        l = summary(l)
        sex_beta = NA
        sex_p = NA
        age_beta = l$coefficients[3,1]
        age_p = l$coefficients[3,4]
      }else{
        l = lm(phenotype ~ variant + sex + age, regression_data)
        l = summary(l)
        sex_beta = l$coefficients[3,1]
        sex_p = l$coefficients[3,4]
        age_beta = l$coefficients[4,1]
        age_p = l$coefficients[4,4]
      }
      
      nm = sum(regression_data$sex == "M")
      nf = sum(regression_data$sex == "F")
      
      variant_beta = l$coefficients[2,1]
      variant_p = l$coefficients[2,4]
      
      hvr_regression_stat = data.frame(
        variant = variant_id,
        phenotype = pheno,
        n_total = n_total,
        n_avail = n_avail,
        nm = nm, 
        nf = nf,
        n0 = n0,
        n1 = n1,
        n2 = n2,
        variant_beta = variant_beta,
        variant_p = variant_p,
        sex_beta = sex_beta,
        sex_p = sex_p,
        age_beta = age_beta,
        age_p = age_p
      )
      
      hvr_regression_stats = rbind.data.frame(
        hvr_regression_stats,
        hvr_regression_stat
      )
    }else{
      cat(paste0("Insufficiant sample size for ", pheno, " to conduct a linear regression...\n"))
    }
    
  }
  
  if(!is.null(dim(sleep_regression_stats))){
    regression_stats = rbind.data.frame(
      mutate(master_regression_stats, table = "master"),
      mutate(hvr_regression_stats, table = "hvr"),
      mutate(sleep_regression_stats, table = "sleep")
    )
  }else{
    regression_stats = rbind.data.frame(
      mutate(master_regression_stats, table = "master"),
      mutate(hvr_regression_stats, table = "hvr")
    )
  }
  
  return(regression_stats)
}

variants = names(select(pheno_master, rs7564509:rs4820434))

variant_stats = vector()
for(v in variants){
  variant_stat = run_regression_by_variant(v)
  variant_stats = rbind.data.frame(variant_stats, variant_stat)
}

variant_stats = variant_stats |>
  arrange(variant_p) |>
  group_by(phenotype, table) |>
  mutate(fdr = p.adjust(variant_p, method = "fdr")) |>
  ungroup() |>
  filter(variant_p <= 0.05)

variant_lookup = variant_lookup %>%
  split(.$variant_id) %>%
  map(function(x){
    data.frame(
      chr = x$chr[1],
      pos = x$pos[1],
      variant = x$variant_id[1],
      gene = paste(x$gene, collapse = ", ")
    )
  }) %>%
  reduce(rbind.data.frame)

variant_stats = variant_stats |>
  left_join(variant_lookup, by = "variant")

write_xlsx(variant_stats, path = "association_stats.xlsx")

