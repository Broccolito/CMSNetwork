library(ggplot2)
library(ggpubr)
library(dplyr)
library(purrr)
library(data.table)
library(scales)
library(patchwork)
library(cowplot)

andean = fread("andean_cms_components.csv")
tibetan = fread("tibetan_cms_components.csv")

color1 = "#788FCE"
color2 = "#E6956F"
threshold_linecolor = "black"

andean = andean %>%
  select(chr, bp, p, cms_score) %>%
  filter(cms_score > 0)

tibetan = tibetan %>%
  select(chr, bp, p, cms_score) %>%
  filter(cms_score > 0)

andean_data = andean
andean_data = andean_data %>%
  mutate(pos = bp) %>%
  select(chr, pos, cms_score) %>%
  mutate(chr = ifelse(chr == "X", "23", chr)) %>%
  mutate(chr = as.numeric(chr)) %>%
  mutate(pos = as.numeric(pos)) %>%
  mutate(cms_score = as.numeric(cms_score))

data_chromosome_position = andean_data %>%
  group_by(chr) %>%
  summarise(chromosome_length = max(pos)) %>%
  mutate(tot = cumsum(as.numeric(chromosome_length)) - chromosome_length) %>%
  select(-chromosome_length)

andean_data_cleaned = left_join(andean_data, data_chromosome_position, by = "chr") %>%
  arrange(chr, pos) %>%
  mutate(bp_cum = pos + tot) %>%
  select(chr, pos, cms_score, tot, bp_cum)

axisdf = andean_data_cleaned %>% 
  group_by(chr) %>% summarize(center = (max(bp_cum) + min(bp_cum))/2)

plt1 = ggplot(andean_data_cleaned, aes(x = bp_cum, y = cms_score)) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c(color1, color2), 22)) +
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
  geom_hline(yintercept = 6, 
             linetype = "dashed", color = threshold_linecolor) +
  xlab(NULL) + 
  ylab(NULL) +
  theme_pubr() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20)
  ); plt1


tibetan_data = tibetan
tibetan_data = tibetan_data %>%
  mutate(pos = bp) %>%
  select(chr, pos, cms_score) %>%
  mutate(chr = ifelse(chr == "X", "23", chr)) %>%
  mutate(chr = as.numeric(chr)) %>%
  mutate(pos = as.numeric(pos)) %>%
  mutate(cms_score = as.numeric(cms_score))

tibetan_data_cleaned = left_join(tibetan_data, data_chromosome_position, by = "chr") %>%
  arrange(chr, pos) %>%
  mutate(bp_cum = pos + tot) %>%
  select(chr, pos, cms_score, tot, bp_cum)

plt2 = ggplot(tibetan_data_cleaned, aes(x = bp_cum, y = -cms_score)) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c(color1, color2), 22)) +
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
  scale_y_continuous(breaks = c(0, -20, -40, -60), labels = c(0, 20, 40, 60)) + 
  geom_hline(yintercept = -6, 
             linetype = "dashed", color = threshold_linecolor) +
  xlab("Chromosome Position") + 
  ylab(NULL) +
  theme_pubr() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 15)
  ); plt2

manh = plt1 / plt2
manh = ggdraw() +
  draw_plot(manh, x = 0.02, y = 0, width = 0.98, height = 1) +
  draw_label("CMS Score", x = 0.02, y = 0.5, angle = 90, vjust = 0.5, hjust = 0.5, size = 15)
manh
ggsave(
  filename = "manhattan_combined.png",
  plot = manh,
  width = 14,        
  height = 6,       
  dpi = 1200,
  bg = "white"
)
