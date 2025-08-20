library(tkoi)
library(dplyr)
library(data.table)
library(writexl)

andean = fread("andean_cms_annotated.csv")
tibetan = fread("tibetan_cms_annotated.csv")

andean = andean |>
  filter(v2g_score > 0.1) |>
  group_by(gene_name) |>
  arrange(desc(cms_score))|>
  slice_max(cms_score, n = 1) |>
  ungroup() |>
  distinct(gene_name, .keep_all = TRUE) |>
  arrange(desc(cms_score))

tibetan = tibetan |>
  filter(v2g_score > 0.1) |>
  group_by(gene_name) |>
  arrange(desc(cms_score))|>
  slice_max(cms_score, n = 1) |>
  ungroup() |>
  distinct(gene_name, .keep_all = TRUE) |>
  arrange(desc(cms_score))

andean_tkoi = andean |>
  select(gene_name, logfc = cms_score, pvalue = p)

tibetan_tkoi = tibetan |>
  select(gene_name, logfc = cms_score, pvalue = p)

set.seed(492357816)
andean_tkoi_result = run_tkoi(
  andean_tkoi,
  subnetwork = tkoi::tkoi_net,
  pvalue_threshold = 0.05,
  logfc_threshold = 0,
  indirect_link_threshold = 3,
  topology_similarity = 0.9,
  n_permutation = 100,
  damping_factor = 0.85,
  maximum_iteration = 500
)

set.seed(492357816)
tibetan_tkoi_result = run_tkoi(
  tibetan_tkoi,
  subnetwork = tkoi::tkoi_net,
  pvalue_threshold = 0.05,
  logfc_threshold = 0,
  indirect_link_threshold = 3,
  topology_similarity = 0.9,
  n_permutation = 100,
  damping_factor = 0.85,
  maximum_iteration = 500
)
 
save(andean_tkoi_result, file = "run_v2g0.1/andean_tkoi_result.rda")
save(tibetan_tkoi_result, file = "run_v2g0.1/tibetan_tkoi_result.rda")

write_xlsx(andean_tkoi_result@network_summary_statistics, 
           path = "run_v2g0.1/andean_tkoi_result.xlsx")
write_xlsx(tibetan_tkoi_result@network_summary_statistics,
           path = "run_v2g0.1/tibetan_tkoi_result.xlsx")






