# Load libraries
library(tidyverse)
library(ggpubr)

# System config
pdf(NULL)

# Load data
basedir = ".."
vaccine_epitope_df = read_tsv(paste0(basedir, "/tables/vaccine_epitope_df.tsv"))
self_similarity_epitopes = read_tsv(paste0(basedir, "/tables/self_similarity_epitopes.tsv"))
epitope_binding_prediction = read.csv(paste0(basedir, "/tables/epitope_binding_prediction_vs_experiment.csv"))

# Prepare data for plotting
## Get Labels
self_similarity_epitopes_labels = self_similarity_epitopes %>%
    group_by(type, med, mea) %>% 
    summarise()

# Plot potential immunogenic epitopes
vaccine_epitope_df %>%
    ggplot(aes(x = vaccine_size, y = epitope_count, group = sample, color = factor(class, levels = rev(c('Missense', 'NOP'))))) +
    geom_line(linewidth = 1.5, alpha = 1) +
    theme_linedraw() +
    scale_color_manual(values = c('#469D88', '#415384')) +
    scale_y_continuous(limits = c(0, 7600), expand = c(0, 0), breaks = c(1000, 3000, 5000, 7000)) +
    scale_x_continuous(limits = c(0, 2000), expand = c(0, 0)) +
    theme(text = element_text(size = 40),
          axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24),
          legend.title = element_blank(),
          legend.position = 'bottom',
          plot.margin = margin(25, 25, 25, 25)) +
    ylab('Possible CD8 Epitopes') +
    xlab('Vaccine size (amino acids)')

ggsave(paste0(basedir, "/figures/epitopes_per_construct.png"), height = 12*0.9, width = 12)

# Plot self similarity boxplots
self_similarity_epitopes %>%
    ggplot(aes(x = reorder(type, med), y = self_sim )) +
    geom_boxplot() +
    theme_linedraw() +
    ylab('Self Similarity') +
    xlab('Epitope Type') +
    stat_compare_means(size = 10, comparisons = list(c('NOP', 'Missense'))) +
    theme(text = element_text(size = 40)) +
    geom_text(data = self_similarity_epitopes_labels, aes(x = type, y = med + 0.02, label = mea), size = 10)

ggsave(paste0(basedir, "/figures/self_sim_boxplot.png"), height = 12*0.8, width = 12)

# Plot epitope binding prediction vs experiment
epitope_binding_prediction %>% 
    ggplot(aes(x = netmhc_filter, y = count, fill = experimental_binder)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste(round(percentage),"%")), position = position_stack(vjust =  0.5)) +
    facet_grid(cols = vars(HLA)) +
    theme_bw() +
    scale_fill_manual(values = c("#E64B35FF", "#00A087FF")) +
    labs(fill = "Experimental Binder") +
    ylab("Epitopes") +
    xlab("netMHCpan EL score filter") +
    theme(text = element_text(size = 16)) +
    theme(legend.position = "bottom")

ggsave(paste0(basedir, "/figures/epitope_binding_prediction.png"), height = 12*0.8, width = 12)