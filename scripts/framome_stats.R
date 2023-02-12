# Load libraries
library(tidyverse)
library(ggsci)
library(patchwork)
library(ggalluvial)

# System config
pdf(NULL)

# Load color pallete 
mypal = pal_npg("nrc", alpha = 1)(10)

# Load data
basedir = ".."
framome_aa = read_tsv(paste0(basedir, "/tables/framome_aa_table.tsv"))
nop_expression_vs_genomic_vaf = read_tsv(paste0(basedir, "/tables/nop_expression_vs_genomic_vaf.tsv"))
sv_flow_df = read_tsv(paste0(basedir, "/tables/sv_flow_table.tsv"))

# Calculate percentage of Genomic Junctions Involved in SV NOP Events
percent_junction = sv_flow_df %>%
    group_by(class, junction) %>%
    summarize() %>%
    count() %>% 
    left_join(
        sv_flow_df %>%
        filter(variants != '') %>%
        group_by(class, junction) %>% 
        summarise() %>% 
        count() %>% 
        rename(a = n)) %>%
    mutate(percent_junction = 100 * a / n) %>%
    mutate(class_mod = case_when(class == "gene_gene" ~ "Gene - Gene",
                                 class == "gene_intergenic" ~ "Gene - Intergenic",
                                 class == "intergenic" ~ "Intergenic",
                                 class == "intragenic" ~ "Intragenic")) %>%
    mutate(class_num = paste0(class_mod, "\n", "n=", a, " (", round(percent_junction, 1), "%)"))

# Order genomic junctions 
ord_percent_junction = c("^Intergenic", "^Intragenic", "^Gene - Intergenic", "^Gene - Gene")
ord_fct_percent_junction = map_chr(ord_percent_junction,
                                   ~percent_junction$class_num[str_detect(percent_junction$class_num, .x)])
percent_junction = percent_junction%>%
mutate(class_num = factor(class_num, levels = ord_fct_percent_junction))%>%
select(-c(n,a,percent_junction, class_mod)) %>%
ungroup()

# Calculate SV NOP Events
sv_nop_events = sv_flow_df %>%
    filter(frame_class != "_") %>%
    group_by(frame_class) %>%
    summarize(count = n()) %>%
    mutate(frame_class_mod = case_when(frame_class == "complex_fusion_gene" ~ "Complex\nGene Fusion",
                                       frame_class == "complex_hidden_frame" ~ "Complex\nHidden NOP",
                                       frame_class == "simple_fusion_gene" ~ "Simple\nGene Fusion",
                                       frame_class == "simple_hidden_frame" ~ "Simple\nHidden NOP")) %>%
    mutate(frame_class_num = paste0(frame_class_mod, "\n", "n=", count)) %>%
    ungroup()

# Order SV NOP Events
ord_sv_nop_events = c("^Complex\nHidden NOP", "^Complex\nGene Fusion", "^Simple\nHidden NOP", "^Simple\nGene Fusion")
ord_fct_sv_nop_events = map_chr(ord_sv_nop_events,
                                ~sv_nop_events$frame_class_num[str_detect(sv_nop_events$frame_class_num, .x)])
sv_nop_events = sv_nop_events %>%
    mutate(frame_class_num = factor(frame_class_num, levels = ord_fct_sv_nop_events)) %>%
    select(-c(count, frame_class_mod))

# Prepare data for plotting
# Add count/precentage
sv_data = left_join(sv_flow_df, percent_junction) %>%
    left_join(sv_nop_events)

# Add multi-NOP; frame_count: 1(Single) or >1(Multiple)
sv_data = sv_data %>%
    mutate(multi_NOP = ifelse(frame_count == 1, 'Single', 'Multiple'))

# Add count to multi-NOP
multi_NOP_num = sv_data %>%
    filter(frame_class != "_") %>%
    group_by(multi_NOP) %>%
    summarize(count = n()) %>%
    mutate(multi_NOP_num = paste0(multi_NOP, "\n", "n=", count)) %>%
    select(-count)
sv_data = left_join(sv_data, multi_NOP_num)

# Add y-axis: Count
sv_data = sv_data %>% mutate(Count = 1)

# Plot framome sizes
textsize = 14
framome_aa %>%
    ggplot(aes(x = reorder(sample, -sum), y = novel_aa, fill = frame_class)) +
    geom_bar(stat = 'identity') +
    ylab("Number of amino acids") +
    facet_grid(. ~ cancer_type, scales = "free", space='free') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = textsize, color="black"),
          legend.position=c(0.5, 0.5), 
          axis.title.x = element_blank(),
          legend.text=element_text(size=textsize, color="black"),
          axis.text.y=element_text(size=textsize, color="black"),
          axis.title.y=element_text(size=textsize, color="black"),
          strip.background = element_rect(color="black", fill="white", linetype="solid"),
          strip.text.x = element_text(size = textsize-2, colour = "black"),
          text = element_text(size = textsize),
          panel.spacing = unit(0.1, "lines")) +
    scale_fill_manual(values = c('Hidden Frame'='#3C5488FF',
                                 'Gene Fusion'='#4DBBD5FF',
                                 'Splice' = '#7E6148FF',
                                 'Indel' = '#E64B35FF',
                                 'Stoploss'='#00A087FF')) +
    labs(fill='NOP category') +
    guides(fill = guide_legend(ncol = 2, title.position='top', title.hjust = 0.5))
       
ggsave(paste0(basedir, "/figures/framome_stats_barplot.png"), height = 12*0.4, width = 12)

# Combined stats plot VAF vs length vs expression
textsize = 22
plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

nop_expression_vs_genomic_vaf_plot <- nop_expression_vs_genomic_vaf %>%
    ggplot(aes(x = tpm, y = vaf, color = frame_class)) +
    geom_point(aes(size = frame_length), position = position_jitter(w = 0.1, h = 0), alpha=0.6) +
    ylab("Purity adjusted VAF") +
    xlab("NOP expression (TPM)") +
    theme_bw() +
    scale_x_log10(labels = plain) +
    theme(legend.position = "bottom", 
          legend.text = element_text(size = textsize, color = "black"),
          legend.title = element_text(size = textsize, color = "black"),
          axis.text.y = element_text(size = textsize, color = "black"),
          axis.title.y = element_text(size = textsize, color = "black"),
          axis.text.x = element_text(size = textsize, color = "black"),
          axis.title.x = element_text(size = textsize, color = "black"),
          strip.background = element_rect(color = "black", fill = "white", linetype = "solid"),
          strip.text.x = element_text(size = textsize, colour = "black")) +
    labs(color = 'NOP category') +
    scale_size(name = "NOP length") +
    guides(colour = guide_legend(override.aes = list(size = 5), ncol = 2, title.position = 'top', title.hjust = 0.5),
           size = guide_legend(override.aes = list(size = 5), ncol = 2, title.position = 'top', title.hjust = 0.5)) +
    scale_color_manual(values = c("Hidden Frame" = "#3C5488FF",
                                  "Gene Fusion" = "#4DBBD5FF",
                                  "Splice" = "#7E6148FF",
                                  "Indel" = "#E64B35FF",
                                  "Stoploss" = "#00A087FF"))

dens1 <- ggplot(nop_expression_vs_genomic_vaf, aes(x = log10(tpm))) + 
  geom_histogram(color = "black", fill = "white") + 
  theme_void()

dens2 <- ggplot(nop_expression_vs_genomic_vaf, aes(x = vaf)) + 
  geom_histogram(color = "black", fill = "white") + 
  theme_void() + 
  coord_flip()

dens1 + plot_spacer() + nop_expression_vs_genomic_vaf_plot + dens2 + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1),
    heights = c(1, 4)) 

ggsave(paste0(basedir, "/figures/framome_stats_point_plot_hist.png"), height = 12, width=12)

# Plot origin of Hidden NOPs and out-of-frame gene fusions
# Theme with no background
theme_nobackground = theme(panel.background = element_blank(),
                           axis.ticks = element_blank(),
                           axis.title.y = element_blank(),
                           axis.text.y = element_blank())

# Color palette
palette_flow = c(mypal[c(10, 9, 3, 2, 5, 6, 8, 4)],
                 "black",
                 "grey")

sv_data %>%
    filter(frame_class != "_") %>%
    select(-c(variants, junction)) %>%
    ggplot(aes(y = Count,
               axis1 = class_num, 
               axis2 = frame_class_num, 
               axis3 = multi_NOP_num)) +
    geom_flow(aes(fill = after_stat(stratum)), alpha = 0.3, width = .5) +
    scale_x_discrete(limits = c("Genomic Junctions\nInvolved in SV NOP Events", 
                                "SV NOP Events", 
                                "Multi-NOP")) +
    geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.7, width = .5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 7.5) +
    scale_fill_manual(values = palette_flow) +
    theme_nobackground +
    theme(legend.position = "none",
          text = element_text(size = 25),
          axis.text.x = element_text(color = "black"))

ggsave(paste0(basedir, "/figures/sv_flow.png"), height = 10, width =17)