library(tkoi)
library(dplyr)
library(data.table)
library(writexl)
library(readxl)
library(igraph)
library(purrr)
library(visNetwork)
library(ggrepel)
library(ggnewscale)

make_overlap_gene_network = function(seed = 492357816){
  
  set.seed(seed)
  gene_overlap = read_excel(path = "network_overlap/Overlapping Gene from Network.xlsx")
  tkoi_net = tkoi::tkoi_net
  
  gene_subg = induced_subgraph(
    tkoi_net,
    vids = V(tkoi_net)[V(tkoi_net)$labels == "['Gene']"]
  )
  
  overlap_gene_nodes = filter(gene_overlap, implicated_in == "Both")$node_id
  neighbor_genes = ego(gene_subg, order = 1, nodes = V(gene_subg)[name %in% overlap_gene_nodes], mode = "all")
  
  neighbor_genes_expanded = map(neighbor_genes, function(x){
    neighbor_nodes = map(x, function(y){
      y[[1]]$name
    }) |>
      unlist()
    names(neighbor_nodes) = NULL
    return(neighbor_nodes)
  })
  names(neighbor_genes_expanded) = names(as.factor(unlist(V(gene_subg)[name %in% overlap_gene_nodes])))
  
  neighbor_genes_expanded = neighbor_genes_expanded |>
    map(function(x){
      gene_overlap |>
        filter(node_id %in% x) |>
        dplyr::select(node_id, name_andean, implicated_in)
    })
  
  neighbor_genes_expanded = do.call(rbind, neighbor_genes_expanded)
  
  
  subg = induced_subgraph(tkoi_net, vids = V(tkoi_net)[name %in% neighbor_genes_expanded$node_id])
  
  # 2) Join attributes (label text, group) onto subgraph vertices
  V(subg)$name_andean  = neighbor_genes_expanded$name_andean[
    match(V(subg)$name, neighbor_genes_expanded$node_id)
  ]
  V(subg)$implicated_in = neighbor_genes_expanded$implicated_in[
    match(V(subg)$name, neighbor_genes_expanded$node_id)
  ]
  
  nodes = data.frame(
    id    = V(subg)$name,
    label = V(subg)$name_andean %||% V(subg)$name, 
    group = V(subg)$implicated_in %||% "Other",
    stringsAsFactors = FALSE
  )
  
  edges = as_data_frame(subg, what = "edges") |>
    dplyr::rename(from = from, to = to)
  
  lay = layout_with_fr(subg, niter = 2000)
  lay = layout.norm(lay, xmin = -300, xmax = 300, ymin = -300, ymax = 300)
  nodes$x = lay[,1]
  nodes$y = lay[,2]
  
  group_cols = c(
    "Andeans-only" = "#9CBEDB",
    "Tibetans-only" = "#48C0AA",
    "Both" = "#EF767A"
  )
  
  group_cols_dark = c(
    "Andeans-only"  = "#1C3B57",   # deep navy
    "Tibetans-only" = "#145247",   # deep teal
    "Both"          = "#66171A"    # deep burgundy
  )
  
  edge_cols = c(
    "Up-regulates"   = "#F6DFD6",
    "Down-regulates" = "#CFCCE3"
  )
  
  edges_xy = edges %>%
    left_join(nodes %>% dplyr::select(id, x, y), by = c("from" = "id")) %>%
    dplyr::rename(x_from = x, y_from = y) %>%
    left_join(nodes %>% dplyr::select(id, x, y), by = c("to" = "id")) %>%
    dplyr::rename(x_to = x, y_to = y) |>
    mutate(edge_kind = ifelse(grepl(pattern = "DOWNREGULATES", edge_type), "Down-regulates", "Up-regulates"))
  
  plt = ggplot() +
    # geom_segment(
    #   data = edges_xy,
    #   aes(x = x_from, y = y_from, 
    #       xend = x_to, yend = y_to),
    #   color = "lightgray",
    #   linewidth = 0.5, alpha = 0.9) +
    geom_segment(
      data = edges_xy,
      aes(x = x_from, y = y_from, 
          xend = x_to, yend = y_to,
          color = edge_kind),
      linewidth = 0.5, alpha = 0.8
    ) +
    scale_color_manual(
      name = "Regulation",
      values = edge_cols
    ) +
    # start a NEW color scale for subsequent layers
    new_scale_color() +
    geom_point(
      data = nodes,
      aes(x = x, y = y, fill = group),
      shape = 21,
      color = "black",
      size = 5
    ) +
    ggrepel::geom_text_repel(
      data = nodes,
      aes(x = x, y = y, label = label, color = group),
      fontface = "italic",
      size = 3, max.overlaps = Inf, 
      box.padding = 0.25, point.padding = 0.2,
      show.legend = FALSE
    ) +
    scale_color_manual(values = group_cols_dark, guide = "none") +
    scale_fill_manual(values = group_cols, name = "Selected in") +  
    coord_equal() +
    theme_void() +
    theme(
      legend.position = "top",
      plot.margin = margin(5, 5, 5, 5)
    )
  
  ggsave(filename = paste0("overlap_gene_network ", seed, ".png"), device = "png", dpi = 1200,
         height = 12, width = 10, bg = "white")
  
}

make_overlap_gene_network(492357816)
