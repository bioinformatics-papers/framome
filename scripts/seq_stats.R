# load libraries
library(tidyverse)
library(cowplot)
library(ggsci)

# System config
pdf(NULL)
textsize = 14

# Palette 
mypal = pal_npg("nrc", alpha = 1)(10)
mypal = mypal[2:10]

# load data
basedir = ".."
somaticvars = read_tsv(paste0(basedir, "/tables/somatic_variants.tsv"))

# Plot panel of sequencing statistics
## SNVs
snv_count = somaticvars %>%
  ggplot(aes(x = reorder(factor(id), -purity), y = snv, fill = cancer_type)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ cancer_type, scales = "free", space = "free") +
  ylab("SNV count") + 
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = textsize, color = "black"),
        axis.title.y = element_text(size = textsize, color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = mypal, name = "My Palette")

## SVs
sv_count = somaticvars %>%
  ggplot(aes(x = reorder(factor(id), -purity), y = sv, fill = cancer_type)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ cancer_type, scales = "free", space = "free") +
  ylab("SV count") +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = textsize, color = "black"),
        axis.title.y = element_text(size = textsize, color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = mypal, name = "My Palette")

## INDELS
indel_count = somaticvars %>%
  ggplot(aes(x = reorder(factor(id), -purity), y = indel, fill = cancer_type)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ cancer_type, scales = "free", space = "free") +
  ylab("Indel count") +
  scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1,
                                   size = textsize,
                                   color = "black"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = textsize, color = "black"),
        axis.title.y = element_text(size = textsize, color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_manual(values = mypal, name = "My Palette")

## PURITY
purity = somaticvars %>%
  ggplot(aes(x = reorder(factor(id), -purity), y = purity, fill = cancer_type)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ cancer_type, scales = "free", space = "free") +
  ylab("Purity (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = textsize, color = "black"),
        axis.title.y = element_text(size = textsize, color = "black"),
        strip.background = element_rect(color = "black",
                                        fill = "white",
                                        linetype = "solid"),
        strip.text.x = element_text(size = textsize, colour = "black")) +
  scale_fill_manual(values = mypal, name = "My Palette")

plot_grid(purity, sv_count, snv_count, indel_count,
          ncol = 1,
          align = "v",
          rel_heights = c(1, 1, 1, 1),
          axis = "lr")

ggsave(paste0(basedir, "/figures/seq_stats_bar_plot.png"), height = 12, width=20)