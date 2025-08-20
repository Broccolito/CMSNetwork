library(dplyr)
library(tkoi)
library(data.table)
library(readxl)
library(writexl)
library(biomaRt)

ensembl = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
andean = fread("andean_cms_annotated.csv")
tibetan = fread("tibetan_cms_annotated.csv")

andean = andean |>
  filter(v2g_score > 0.2) |>
  group_by(gene_name) |>
  arrange(desc(cms_score))|>
  slice_max(cms_score, n = 1) |>
  ungroup() |>
  distinct(gene_name, .keep_all = TRUE) |>
  arrange(desc(cms_score))

tibetan = tibetan |>
  filter(v2g_score > 0.2) |>
  group_by(gene_name) |>
  arrange(desc(cms_score))|>
  slice_max(cms_score, n = 1) |>
  ungroup() |>
  distinct(gene_name, .keep_all = TRUE) |>
  arrange(desc(cms_score))

gene_map = getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(c(andean$gene_name, tibetan$gene_name)),
  mart = ensembl
)

overlap = inner_join(andean, tibetan, by = "gene_name", suffix = c("_andean", "_tibetan")) |>
  left_join(gene_map, by = c("gene_name" = "ensembl_gene_id"))

overlap = overlap |>
  dplyr::select(gene_name, hgnc_symbol, chr = chr_andean, everything()) |>
  dplyr::select(-chr_tibetan)

write_xlsx(overlap, path = "run_v2g0.2/tibetan_andean_overlap.xlsx")



