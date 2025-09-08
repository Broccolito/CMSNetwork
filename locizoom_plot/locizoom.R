locizoom = function(gene_name, rs_num,
                    snp_chr, snp_loc,
                    threshold, file_name,
                    plot_ratio = 0.3,
                    plot_height = 8, 
                    plot_width = 7,
                    text_size = 15,
                    text_shift = 0.4,
                    text_horizontal_shift = 0,
                    include_genetrack = FALSE,
                    include_tibetan = TRUE,
                    andean_color = "#D39200",
                    tibetan_color = "#7CAE00",
                    save_plot = TRUE,
                    return_markers = TRUE){
  
  lower = snp_loc - threshold
  upper = snp_loc + threshold
  
  cms_andeans = cms_hg19$andean_cms %>%
    filter(chr == snp_chr) %>%
    filter(pos >= lower) %>%
    filter(pos <= upper) %>%
    mutate(pos = pos / 1e6)
  plot_xrange = range(cms_andeans$pos)
  
  cms_tibetans = cms_hg19$tibetan_cms %>%
    filter(cms_score > 0) %>%
    filter(chr == snp_chr) %>%
    filter(pos >= lower) %>%
    filter(pos <= upper) %>%
    mutate(pos = pos / 1e6)
  
  if(dim(cms_tibetans)[1]==0){
    include_tibetan = FALSE
    cat("There is no Tibetan CMS marker within the region...\n")
  }
  
  recomb = cms_hg19$recombination_rate %>%
    filter(chr == paste0("chr",snp_chr)) %>%
    filter(pos >= lower) %>%
    filter(pos <= upper) %>%
    mutate(pos = pos / 1e6)
  
  genes = cms_hg19$gene_annotation %>%
    filter(Chrom == snp_chr) %>%
    filter(End >= lower) %>%
    filter(Start <= upper) %>%
    mutate(Start = Start / 1e6) %>%
    mutate(End = End / 1e6)
  genes = genes %>%
    mutate(index = seq(1:dim(genes)[1])) %>%
    mutate(Start = ifelse(Start<=plot_xrange[1],plot_xrange[1],Start)) %>%
    mutate(End = ifelse(End>=plot_xrange[2],plot_xrange[2],End))
  
  is_cms_marker = (snp_loc %in% cms_hg19$andean_cms$pos)
  snp_location = tibble(x = snp_loc, y = ifelse((snp_loc %in% cms_hg19$andean_cms$pos),
                                                filter(cms_hg19$andean_cms,pos==snp_loc)$cms_score,
                                                6),
                        text = paste0(gene_name," ",rs_num)) %>%
    mutate(x = x / 1e6)
  
  if(include_tibetan){
    multiplier = 100 / max(c(cms_andeans$cms_score, cms_tibetans$cms_score))
  }else{
    multiplier = 100 / max(c(cms_andeans$cms_score))
  }
  
  recomb$rate = recomb$rate / multiplier
  
  
  if(include_tibetan){
    plt1 = ggplot() + 
      geom_line(data = recomb, aes(x = pos, y = rate),
                color = "#0080FE", linewidth = 0.5) + 
      geom_point(data = cms_andeans, 
                 aes(x = pos,y = ifelse(cms_score==0,NA,cms_score)), 
                 size = 3, color = "black", fill = andean_color,
                 shape = 21) +
      geom_point(data = cms_tibetans, 
                 aes(x = pos,y = ifelse(cms_score==0,NA,cms_score)), 
                 size = 3, color = "black", fill = tibetan_color,
                 shape = 21) +
      geom_hline(yintercept = 6, color = "black",linetype = 2) + 
      # geom_point(data = snp_location,
      #            aes(x = x, y = y),
      #            color = "black", fill = "red",
      #            size = 4, shape = ifelse(is_cms_marker,21,22)) + 
      geom_text(data = snp_location,
                aes(x = x+text_horizontal_shift, y = y+text_shift, label = text)) + 
      scale_y_continuous(
        name = TeX("$CMS_{BF}$"),
        sec.axis = sec_axis(transform = ~.*multiplier, 
                            name = "Recombination Rate (cM/Mb)")
        
      ) + 
      xlab("Chromosome Position (Mb)") + 
      theme_bw() + 
      theme(text = element_text(size = text_size),
            axis.title.y.right = element_text(color = "#0080FE"))
  }else{
    plt1 = ggplot() + 
      geom_line(data = recomb, aes(x = pos, y = rate),
                color = "#0080FE", linewidth = 0.5) + 
      geom_point(data = cms_andeans, 
                 aes(x = pos,y = ifelse(cms_score==0,NA,cms_score)), 
                 size = 3, color = "black", fill = andean_color,
                 shape = 21) +
      geom_hline(yintercept = 6, color = "black",linetype = 2) + 
      # geom_point(data = snp_location,
      #            aes(x = x, y = y),
      #            color = "black", fill = "red",
      #            size = 4, shape = ifelse(is_cms_marker,21,22)) + 
      geom_text(data = snp_location,
                aes(x = x+text_horizontal_shift, y = y+text_shift, label = text)) + 
      scale_y_continuous(
        name = TeX("$CMS_{BF}$"),
        sec.axis = sec_axis(transform = ~.*multiplier, 
                            name = "Recombination Rate (cM/Mb)")
      ) + 
      xlab("Chromosome Position (Mb)") + 
      theme_bw() + 
      theme(text = element_text(size = text_size),
            axis.title.y.right = element_text(color = "#0080FE"))
  }
  
  
  plt2 = ggplot() + 
    geom_errorbarh(data = genes,
                   aes(y = index,xmax = End, xmin = Start,
                       height = 0.1, color = Coding,
                       group = Gene),
                   linewidth = 1) +
    geom_text(data = genes,
              aes(x = (Start+End)/2, y = index + 0.5,label = Gene)) + 
    geom_vline(xintercept = plot_xrange[1]) + 
    geom_vline(xintercept = plot_xrange[2]) + 
    xlim(plot_xrange[1],plot_xrange[2]) +
    xlab("") + 
    ylab("") + 
    labs(color = "") + 
    theme_void() + 
    theme(text = element_text(size = text_size),
          axis.text.y = element_text(color = "white"),
          axis.ticks = element_blank(),
          legend.position = "bottom"
    ) 
  
  if(include_genetrack){
    suppressWarnings({
      plt = ggarrange(plt1,plt2, ncol = 1,heights = c(1, plot_ratio))
    })
  }else{
    plt = plt1
  }
  
  if(save_plot){
    graphics.off()
    for(i in 1:2){
      suppressWarnings({
        try({
          ggsave(filename = paste0("saved_figures/",file_name), 
                 plot = plt, device = "png",
                 dpi = 1200, height = plot_height,width = plot_width)
          
        },silent = TRUE)
      })
    }
  }

  if(return_markers){
    cms_andeans = cms_andeans %>%
      filter(cms_score > 0) %>%
      mutate(ancestry = "Andeans") %>%
      select(chr, pos, cms_score, ancestry)
    cms_tibetans = cms_tibetans %>%
      filter(cms_score > 0) %>%
      mutate(ancestry = "Tibetans") %>%
      select(chr, pos, cms_score, ancestry)
    
    if(dim(cms_tibetans)[1]==0){
      all_markers = cms_andeans
    }else{
      if(dim(cms_andeans)[1]==0){
        all_markers = cms_tibetans
      }else{
        all_markers = rbind.data.frame(cms_andeans, cms_tibetans)
      }
    }
    
    if(dim(all_markers)[1]!=0){
      all_markers = all_markers %>%
        mutate(cms_marker = paste("chr",chr,":",pos*1e6, sep = "")) %>%
        mutate(alias = cms_marker)
    }
    
    return(all_markers)
    
  }else{
    
    return(plt)
    
  }
  
}



