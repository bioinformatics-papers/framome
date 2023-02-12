# Load libraries
library(tidyverse)
library(ggsci)

# System config
pdf(NULL)

# Load data
basedir = ".."
riboseq_coverage = read_tsv(paste0(basedir, "/tables/hidden_frame_riboseq_coverage.tsv"))

# Plot Ribo-seq segnaling for hidden NOPs
riboseq_coverage %>% 
    group_by(frame_peptide) %>% 
    mutate(total_p_site_coverage = sum(riboseq_reads)) %>% 
    ggplot(aes(x = reorder(frame_peptide, -total_p_site_coverage), y = riboseq_reads, fill = reading_frame)) +
    geom_bar(stat = "identity") +
    theme_linedraw() +
    facet_grid(cols = vars(sample), scales = "free_x", space = "free") +
    scale_fill_npg() +
    theme(legend.position = c(0.6, 0.6)) +
    labs(fill='Phase') +
    theme(text = element_text(size = 18)) +
    theme(axis.text.x = element_blank()) +
    xlab("Hidden Frame") +
    ylab("RiboSeq\nFragments") +
    theme(strip.text = element_text(colour = "black", size = 18)) +
    theme(strip.background = element_rect(fill = NA))

ggsave(paste0(basedir, "/figures/riboseq_signaling_hidden_NOPs.png"), height = 12*0.3, width = 12)