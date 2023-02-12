# Load libraries
library(tidyverse)
library(ggsci)
library(scales)

# System config
pdf(NULL)

# Load data
basedir = ".."
binding_data = read_tsv(paste0(basedir, "/tables/summarized_lung_binding.tsv"))
jitter_binding_data = read_tsv(paste0(basedir, "/tables/summarized_lung_binding_jitter.tsv"))

# Prepare data to plot
## Get peptide order
peptide_order = (binding_data %>%
    group_by(peptide) %>%
    summarize(r = max(r)) %>%
    arrange(-r))$peptide

## Get HLA order
hla_order = rev((binding_data %>%
    group_by(hla) %>%
    summarize %>%
    arrange(hla))$hla)

# Plot binding assays
binding_data %>%
    ggplot(aes(x = factor(peptide, levels = peptide_order), y = avg, fill = factor(hla, levels = hla_order))) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw() +
    facet_grid(cols = vars(sample), scales = 'free_x', space = 'free') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_npg() +
    geom_jitter(aes(y = value), alpha = 1, height = 0, width = 0.1, data = jitter_binding_data) +
    geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd)) +
    theme(legend.title = element_blank(),
          legend.position = 'right',
          axis.title.x = element_blank(),
          text= element_text(size = 14)) +
    ylab('Percentage binding relative\nto positive control') +
    scale_y_continuous(breaks = breaks_width(20))

ggsave(paste0(basedir, "/figures/lung_binders.png"), height = 12*0.3, width = 12)