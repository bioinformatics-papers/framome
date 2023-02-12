# Load libraries
library(tidyverse)
library(ggsci)
library(scales)
library(cowplot)
set.seed(7272)

# System config
pdf(NULL)

# load nature palette 
mypal = pal_npg("nrc", alpha = 1)(10)
mypal = mypal[2:10]

# Load lrRNA read length data
basedir = ".."
lrRNA_read_length = read_tsv(paste0(basedir, "/tables/lrRNA_read_length.tsv"))
lrRNA_number_reads = read_tsv(paste0(basedir, "/tables/lrRNA_number_reads.tsv"))
lrRNA_full_length = read_tsv(paste0(basedir, "/tables/full_length_fraction_lrRNA.csv"))

# Read length boxplot
read_length_boxplot = lrRNA_read_length %>%
    ggplot(aes(reorder(factor(Alias), -Purity), Lengths, fill = Cancer_Type)) +
    geom_boxplot() +
    geom_point(aes(reorder(factor(Alias), -Purity), N50),
              size = 4, shape = 24, color = "black", fill = "red") +
    coord_cartesian(ylim = c(0,2500)) +
    scale_y_continuous(labels = comma) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 25, color = "black"),
          axis.title.x = element_text(size = 30, color = "black", margin = margin(t = 20)),
          axis.text.y = element_text(size = 25, color = "black"),
          axis.title.y = element_text(size = 30, color = "black", margin = margin(r = 20)),
          strip.background = element_rect(color = "black", fill = "white", linetype = "solid"),
          strip.text.x = element_text(size = 28, color = "black"),
          plot.title = element_text(size = 35)) +
    facet_grid(. ~ Cancer_Type, scales = "free", space = 'free') +
    scale_fill_npg() +
    xlab("Samples") +
    ylab("Read length (bp)") +
    ggtitle("Read lengths (Random 100K) ordered by cancer type and purity(+ to -)") +
    scale_fill_manual(values = mypal, name = "My Palette")

# Number of reads barplot
number_reads = lrRNA_number_reads %>%
    ggplot(aes(reorder(factor(Alias), -Purity), Reads, fill = Cancer_Type)) +
    geom_col() +
    scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6), breaks = seq(0, 140000000, by = 20000000)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 25, color = "black"),
          axis.title.x = element_text(size = 30, color = "black", margin = margin(t = 20)),
          axis.text.y = element_text(size = 25, color = "black"),
          axis.title.y = element_text(size = 30, color = "black", margin = margin(r = 20)),
          strip.background = element_rect(color = "black", fill = "white", linetype = "solid"),
          strip.text.x = element_text(size = 28, color = "black"),
          plot.title = element_text(size = 35)) +
    facet_grid(. ~ Cancer_Type, scales = "free", space = 'free') +
    scale_fill_npg() +
    xlab("Samples") +
    ylab("Number of reads (millions)") +
    ggtitle("Number of reads ordered by cancer type and purity(+ to -)") +
    scale_fill_manual(values = mypal, name = "My Palette")

# Number of full-length reads barplot
full_length = lrRNA_full_length %>%
    ggplot(aes(reorder(factor(Alias), -Purity), alias_fraction_full_length, fill = Cancer_Type)) +
    geom_col() +
    scale_y_continuous(labels = percent_format(), limits = c(0,1)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 25, color = "black"),
          axis.title.x = element_text(size = 30, color = "black",margin = margin(t = 20)),
          axis.text.y = element_text(size = 25, color = "black"),
          axis.title.y = element_text(size = 30, color = "black", margin = margin(r = 20)),
          strip.background = element_rect(color = "black", fill = "white", linetype = "solid"),
          strip.text.x = element_text(size = 28, color = "black"),
          plot.title = element_text(size = 35)) +
    facet_grid(. ~ Cancer_Type, scales = "free", space = 'free') +
    scale_fill_npg() +
    xlab("Samples") +
    ylab("Full-length reads (%)") +
    ggtitle("Number of full-length reads in protein coding genes ordered by cancer type and purity(+ to -)") +
    scale_fill_manual(values = mypal, name = "My Palette")

# Remove elements from individual plots
grid_theme = theme(plot.title = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
no_Xheader = theme(strip.text.x = element_blank())

read_length_boxplot = read_length_boxplot + grid_theme
number_reads = number_reads + grid_theme + no_Xheader
full_length = full_length + theme(plot.title = element_blank()) + no_Xheader

# Create a panel with read length, number of reads and full-length reads
plot_grid(read_length_boxplot, number_reads, full_length, ncol = 1, 
          align = "v", 
          rel_heights=c(1,1,1),
          axis = 'lr')

ggsave(paste0(basedir, "/figures/lrRNA_seq_stats_panel.png"), height = 20, width=35)