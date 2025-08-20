library(dplyr)
library(arrow)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(ggpubr)
library(RNOmni)
library(biomaRt)
library(readxl)
library(writexl)

mart = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes_map = getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(overlap$gene_name),
  mart = mart
)

v2g = fread("v2g.csv.gz")

tibetan = fread("tibetan_cms_components.csv")
andean = fread("andean_cms_components.csv")

run_liftover = function(dataset){
  chain = import.chain("hg19ToHg38.over.chain")
  dataset_liftover = GRanges(seqnames = paste0("chr", dataset$chr),
                             ranges = IRanges(start = dataset$bp, end = dataset$bp)) |>
    liftOver(chain = chain)
  
  dataset = dataset[lengths(dataset_liftover) == 1]
  dataset_liftover = unlist(dataset_liftover[lengths(dataset_liftover) == 1])
  
  dataset$chr_hg38 = as.character(seqnames(dataset_liftover))
  dataset$bp_hg38  = start(dataset_liftover)
  
  return(dataset)
}

tibetan = tibetan |>
  mutate(chr = as.character(chr))

andean = andean |>
  mutate(chr = as.character(chr))

tibetan = run_liftover(tibetan)
andean = run_liftover(andean)

tibetan_annotated = v2g |>
  inner_join(tibetan, by = c("chr" = "chr", "position" = "bp_hg38"))

andean_annotated = v2g |>
  inner_join(andean, by = c("chr" = "chr", "position" = "bp_hg38"))

tibetan_annotated = tibetan_annotated |>
  mutate(chr_hg38 = gsub(pattern = "chr", replacement = "", chr_hg38)) |>
  filter(chr == chr_hg38) |>
  dplyr::select(chr, bp, position, ref_allele, alt_allele,
         gene_id, overall_score, 
         cms, cms_score, p, xpehh, ihs, dihh, fst, daf)
names(tibetan_annotated) = c("chr", "pos_hg37", "pos_hg38", "ref", "alt", 
                             "gene_name", "v2g_score",
                             "cms", "cms_score", "p", "xpehh", "ihs", "dihh", "fst", "daf")

andean_annotated = andean_annotated |>
  mutate(chr_hg38 = gsub(pattern = "chr", replacement = "", chr_hg38)) |>
  filter(chr == chr_hg38) |>
  dplyr::select(chr, bp, position, ref_allele, alt_allele,
         gene_id, overall_score, 
         cms, cms_score, p, xpehh, ihs, dihh, fst, daf)
names(andean_annotated) = c("chr", "pos_hg37", "pos_hg38", "ref", "alt", 
                            "gene_name", "v2g_score",
                            "cms", "cms_score", "p", "xpehh", "ihs", "dihh", "fst", "daf")

andean_annotated = andean_annotated |>
  filter(v2g_score > 0.2) |>
  group_by(gene_name) |>
  arrange(desc(cms_score))|>
  slice_max(cms_score, n = 1) |>
  ungroup() |>
  distinct(gene_name, .keep_all = TRUE) |>
  arrange(desc(cms_score))

tibetan_annotated = tibetan_annotated |>
  filter(v2g_score > 0.2) |>
  group_by(gene_name) |>
  arrange(desc(cms_score))|>
  slice_max(cms_score, n = 1) |>
  ungroup() |>
  distinct(gene_name, .keep_all = TRUE) |>
  arrange(desc(cms_score))

overlap = inner_join(andean_annotated, tibetan_annotated,
                     by = "gene_name", 
                     suffix = c("_andean", "_tibetan")) |>
  mutate(transformed_beta_tibetan = RankNorm(cms_score_tibetan)) |>
  mutate(transformed_beta_andean = RankNorm(cms_score_andean))

overlap = overlap |>
  left_join(genes_map, by = c("gene_name" = "ensembl_gene_id")) |>
  mutate(gene_symbol = hgnc_symbol) |>
  dplyr::select(-hgnc_symbol) |>
  dplyr::select(chr_andean:gene_name, gene_symbol, everything())

write_xlsx(overlap, path = "all_tibetan_andean_overlap.xlsx")


overlap = overlap |>
  mutate(point_color = ifelse(cms_score_tibetan >= 6 & cms_score_andean >= 6, "#EF767A", ifelse(
    cms_score_andean >= 6, "#9CBEDB", ifelse(
    cms_score_tibetan >= 6, "#48C0AA", "gray"
  )
))) 

overlap_significant = filter(overlap, cms_score_andean >= 6 & cms_score_tibetan >= 6)

write_xlsx(overlap_significant, path = "tibetan_andean_overlap_significant.xlsx")

lm_stats = lm(formula = transformed_beta_andean ~ transformed_beta_tibetan, overlap) |>
  summary()
lm_p = as.character(signif(lm_stats$coefficients[2,4], digits = 3))
lm_p = ifelse(lm_p == "0", "1e-200", lm_p)
lm_slope = as.character(round(lm_stats$coefficients[2,1], 2))
rsq = as.character(round(lm_stats$r.squared, 3))

plt_raw = ggplot(data = overlap, aes(x = transformed_beta_tibetan,
                                 y = transformed_beta_andean)) + 
  geom_point(color = overlap$point_color) +
  geom_point(data = overlap_significant,
             size = 2, color = "black", shape = 21, fill = overlap_significant$point_color) + 
  geom_smooth(method = "lm", formula = "y~x", color = "black") +
  xlab(expression("Transformed " * CMS[BF] * " in Tibetans")) +
  ylab(expression("Transformed " * CMS[BF] * " in Andeans")) +
  ggtitle(label = "Whole-Genome Selection Scan",
          subtitle = bquote("Slope = " * .(lm_slope) * 
                              "; P.Value < " * .(lm_p) * 
                              "; " * R^2 * " = " * .(rsq))) +
  theme_bw() + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))

plot_scatter = function(modality){
  tibetan_genes = read_excel(path = "run_v2g0.2/tibetan_tkoi_result.xlsx", sheet = modality)
  andean_genes = read_excel(path = "run_v2g0.2/andean_tkoi_result.xlsx", sheet = modality)
  
  overlap = inner_join(tibetan_genes, andean_genes, 
                       suffix = c("_tibetan", "_andean"),
                       by = "node_id") |>
    mutate(fdr_tibetan = ifelse(fdr_tibetan == 0, 1e-200, fdr_tibetan)) |>
    mutate(fdr_andean = ifelse(fdr_andean == 0, 1e-200, fdr_andean)) |>
    filter(!is.na(fdr_andean)) |>
    filter(!is.na(fdr_tibetan)) |>
    mutate(point_label = ifelse(fdr_andean <= 0.05 & fdr_tibetan <= 0.05, name_andean, NA)) |>
    mutate(point_color = ifelse(fdr_andean <= 0.05 & fdr_tibetan <= 0.05, "#EF767A", ifelse(
      fdr_andean <= 0.05, "#9CBEDB", ifelse(
        fdr_tibetan <= 0.05, "#48C0AA", "gray"
      )
    ))) |>
    mutate(transformed_beta_tibetan = RankNorm(beta_tibetan)) |>
    mutate(transformed_beta_andean = RankNorm(beta_andean))
  
  overlap_significant = filter(overlap, fdr_andean <= 0.05 & fdr_tibetan <= 0.05)
  
  lm_stats = lm(formula = transformed_beta_andean ~ transformed_beta_tibetan, overlap) |>
    summary()
  lm_p = as.character(signif(lm_stats$coefficients[2,4], digits = 3))
  lm_slope = as.character(round(lm_stats$coefficients[2,1], 2))
  rsq = as.character(round(lm_stats$r.squared, 3))
  
  plt = ggplot(data = overlap, aes(x = transformed_beta_tibetan,
                                   y = transformed_beta_andean)) + 
    geom_point(color = overlap$point_color) +
    geom_point(data = overlap_significant,
               size = 2, color = "black", shape = 21, fill = overlap_significant$point_color) + 
    geom_smooth(method = "lm", formula = "y~x", color = "black") +
    xlab(expression("Transformed " * beta * " in Tibetans")) +
    ylab(expression("Transformed " * beta * " in Andeans")) +
    ggtitle(label = modality,
            subtitle = bquote("Slope = " * .(lm_slope) * 
                                "; P.Value = " * .(lm_p) * 
                                "; " * R^2 * " = " * .(rsq))) +
    theme_bw() + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
  
  return(plt)
}

plt_network = plot_scatter("Gene")

panel = (plt_raw | plt_network) +
  plot_annotation(tag_levels = 'A')

ggsave("gene_level_convergence.png", plot = panel, 
       width = 8, height = 4, dpi = 1200)


