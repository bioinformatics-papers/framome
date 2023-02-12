# Load libraries
library(tidyverse)

# System config
pdf(NULL)
textsize = 14

# Load data
basedir = ".."
all_tpm_table = read_tsv(paste0(basedir, "/tables/all_tpm_table.tsv"))
nop_vs_missense = read_csv(paste0(basedir, "/tables/nop_vs_missense.csv"))
nop_purity_vaf_adjusted = read_csv(paste0(basedir, "/tables/nop_purity_vaf_adjusted.csv"))

# Plot framome percentiles 
all_tpm_table %>%
  mutate(frame_class = factor(frame_class, levels = c('Stoploss', 'Splice', 'Indel', 'Gene Fusion', 'Hidden NOP'))) %>%
  ggplot(aes(x = reorder(sample, -median), y = 100 * sample_percentile)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0, aes(color = frame_class)) +
  ylab("Protein Coding\nExpression Percentile") +
  facet_grid(. ~ cancer_type, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = textsize, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = textsize, color = "black"),
        axis.text.y = element_text(size = textsize, color = "black"),
        axis.title.y = element_text(size = textsize, color = "black"),
        strip.background = element_rect(color = "black", fill = "white", linetype = "solid"),
        strip.text.x = element_text(size = textsize - 2, color = "black"),
        text = element_text(size = textsize),
        legend.position = "bottom",
        panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = c('Hidden NOP' = '#3C5488FF',
                               'Gene Fusion' = '#4DBBD5FF',
                               'Splice' = '#7E6148FF',
                               'Indel' = '#E64B35FF',
                               'Stoploss' = '#00A087FF')) +
  labs(color = "NOP category") +
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size = 5), ncol = 5, title.position = "top", title.hjust = 0.5))

ggsave(paste0(basedir, "/figures/frame_expression_percentiles.png"), height = 12*0.4, width = 12)

# Plot frame RNA expression
nop_vs_missense %>%
    ggplot(aes(x = factor(frame_class, levels = c("Normal\nGene", 
                                                  "Missense",
                                                  "Gene\nFusion",
                                                  "Hidden\nNOP",
                                                  "Indel",
                                                  "Stoploss",
                                                  "Splice")),
               y = tpm)) +
    geom_boxplot() +
    theme_bw() +
    coord_cartesian(ylim = c(0, 17)) +
    theme(text = element_text(size = textsize),
          axis.title.x = element_blank()) +
    ylab('RNA Expression (TPM)')

ggsave(paste0(basedir, "/figures/frame_expression_boxplots.png"), height = 6, width=7)

# Frame VAF boxplots
nop_purity_vaf_adjusted %>%
    ggplot(aes(x = factor(frame_class, levels = c("Stoploss",
                                                  "Splice",
                                                  "Missense",
                                                  "Indel",
                                                  "Gene\nFusion",
                                                  "Hidden\nNOP")),
               y = vaf)) +
    geom_boxplot() +
    theme_bw() +
    coord_cartesian(ylim = c(0, 1)) +
    theme(text = element_text(size = textsize),
          axis.title.x = element_blank()) +
    ylab('Purity adjusted VAF')

ggsave(paste0(basedir, "/figures/frame_vaf_boxplots.png"), height = 6, width=7)