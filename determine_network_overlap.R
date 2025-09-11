library(dplyr)
library(tkoi)
library(data.table)
library(readxl)
library(writexl)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(RNOmni)
library(patchwork)

calc_stats = function(modality){
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
  lm_p = signif(lm_stats$coefficients[2,4], digits = 3)
  lm_slope = round(lm_stats$coefficients[2,1], 2)
  rsq = round(lm_stats$r.squared, 3)
  
  stats = data.frame(
    modality = modality,
    lm_p = lm_p,
    lm_slope = lm_slope,
    rsq = rsq
  )
  return(stats)
}

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

plt3 = plot_scatter("Anatomy")
plt1 = plot_scatter("CellType")
plt6 = plot_scatter("Disease")
plt7 = plot_scatter("BiologicalProcess")
plt8 = plot_scatter("MolecularFunction")
plt9 = plot_scatter("CellularComponent")
plt5 = plot_scatter("ProteinDomain")
plt4 = plot_scatter("Pathway")
plt2 = plot_scatter("Reaction")

strip_x = theme(axis.title.x = element_blank())
strip_y = theme(axis.title.y = element_blank())

plt1 = plt1 + strip_x + strip_y
plt2 = plt2 + strip_x + strip_y
plt3 = plt3 + strip_x + strip_y
plt4 = plt4 + strip_x
plt5 = plt5 + strip_x + strip_y
plt6 = plt6 + strip_x + strip_y
plt7 = plt7 + strip_x + strip_y
plt8 = plt8 + strip_y
plt9 = plt9 + strip_x + strip_y

panel = (plt1 | plt2 | plt3) /
  (plt4 | plt5 | plt6) /
  (plt7 | plt8 | plt9)

ggsave("network_convergence.png", plot = panel, 
       width = 12, height = 8, dpi = 1200)


overlap_stats = rbind.data.frame(
  calc_stats("Anatomy"),
  calc_stats("CellType"),
  calc_stats("Disease"),
  calc_stats("BiologicalProcess"),
  calc_stats("MolecularFunction"),
  calc_stats("CellularComponent"),
  calc_stats("ProteinDomain"),
  calc_stats("Pathway"),
  calc_stats("Reaction")
) |>
  arrange(desc(lm_slope))


overlap_stats_ordered = overlap_stats[order(overlap_stats$lm_slope, decreasing = TRUE), ]
plot_data = melt(overlap_stats_ordered[, c("modality", "lm_slope", "rsq")], 
                 id.vars = "modality")

plot_data$pos = ifelse(plot_data$variable == "lm_slope", -0.2, 0.2)
scale_factor = max(overlap_stats_ordered$lm_slope) / max(overlap_stats_ordered$rsq)

barplot_figure = ggplot(plot_data, aes(x = factor(modality, levels = overlap_stats_ordered$modality))) +
  geom_col(data = subset(plot_data, variable == "lm_slope"), color = "black",
           aes(y = value), fill = "steelblue", width = 0.4, position = position_nudge(x = -0.2)) +
  geom_col(data = subset(plot_data, variable == "rsq"), color = "black",
           aes(y = value * scale_factor), fill = "coral", width = 0.4, position = position_nudge(x = 0.2)) +
  scale_y_continuous(
    name = expression(paste("Regression Slope ")),
    sec.axis = sec_axis(~ . / scale_factor, name = expression(paste("Regression ", R^2)))
  ) +
  labs(title = NULL, x = NULL) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.left = element_text(color = "steelblue"),
    axis.title.y.right = element_text(color = "coral")
  )

ggsave("regression_stats_barplot.png", plot = barplot_figure, 
       width = 12, height = 3.5, dpi = 1200)

export_overlap = function(modality = "Gene"){
  
  tibetan_genes = read_excel(path = "run_v2g0.2/tibetan_tkoi_result.xlsx", sheet = modality)
  andean_genes = read_excel(path = "run_v2g0.2/andean_tkoi_result.xlsx", sheet = modality)
  
  overlap = inner_join(tibetan_genes, andean_genes, 
                       suffix = c("_tibetan", "_andean"),
                       by = "node_id") |>
    mutate(fdr_tibetan = ifelse(fdr_tibetan == 0, 1e-200, fdr_tibetan)) |>
    mutate(fdr_andean = ifelse(fdr_andean == 0, 1e-200, fdr_andean)) |>
    filter(!is.na(fdr_andean)) |>
    filter(!is.na(fdr_tibetan)) |>
    mutate(implicated_in = ifelse(fdr_andean <= 0.05 & fdr_tibetan <= 0.05, "Both", ifelse(
      fdr_andean <= 0.05, "Andeans-only", ifelse(
        fdr_tibetan <= 0.05, "Tibetans-only", "Neither"
      )
    ))) |>
    mutate(transformed_beta_tibetan = RankNorm(beta_tibetan)) |>
    mutate(transformed_beta_andean = RankNorm(beta_andean))
  
  overlap$implicated_in = factor(overlap$implicated_in, 
                                 levels = c("Both", "Andeans-only", "Tibetans-only", "Neither"))
  
  overlap_significant = filter(overlap, fdr_andean <= 0.05 | fdr_tibetan <= 0.05) |>
    arrange(implicated_in)
  
  write_xlsx(overlap_significant, path = file.path("network_overlap", 
                                                   paste0("Overlapping ", modality, " from Network.xlsx")))
  
}


export_overlap("Anatomy")
export_overlap("CellType")
export_overlap("Disease")
export_overlap("BiologicalProcess")
export_overlap("MolecularFunction")
export_overlap("CellularComponent")
export_overlap("ProteinDomain")
export_overlap("ProteinFamily")
export_overlap("Protein")
export_overlap("Complex")
export_overlap("Pathway")
export_overlap("Reaction")
export_overlap("Gene")


