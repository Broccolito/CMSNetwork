library(dplyr)
library(arrow)
library(data.table)
library(rtracklayer)
library(GenomicRanges)

v2g = fread("v2g.csv")

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
  filter(cms_score > 6) |>
  mutate(chr = as.character(chr))

andean = andean |>
  filter(cms_score > 6) |>
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
  select(chr, bp, position, ref_allele, alt_allele,
         gene_id, overall_score, 
         cms, cms_score, p, xpehh, ihs, dihh, fst, daf)
names(tibetan_annotated) = c("chr", "pos_hg37", "pos_hg38", "ref", "alt", 
                             "gene_name", "v2g_score",
                             "cms", "cms_score", "p", "xpehh", "ihs", "dihh", "fst", "daf")

andean_annotated = andean_annotated |>
  mutate(chr_hg38 = gsub(pattern = "chr", replacement = "", chr_hg38)) |>
  filter(chr == chr_hg38) |>
  select(chr, bp, position, ref_allele, alt_allele,
         gene_id, overall_score, 
         cms, cms_score, p, xpehh, ihs, dihh, fst, daf)
names(andean_annotated) = c("chr", "pos_hg37", "pos_hg38", "ref", "alt", 
                             "gene_name", "v2g_score",
                             "cms", "cms_score", "p", "xpehh", "ihs", "dihh", "fst", "daf")

fwrite(tibetan_annotated, file = "tibetan_cms_annotated.csv")
fwrite(andean_annotated, file = "andean_cms_annotated.csv")

