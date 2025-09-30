# relative abundance microViz plots
# 16S only
all_16s_phyloseq_rare_merged <- merge_treatments(all_16s_phyloseq_rare, c("Country", "Depth_cm", "CNP_Scenario"))

all_16s_phyloseq_rare_merged_phyla <- all_16s_phyloseq_rare_merged %>%
  ps_select(Country, Depth_cm, CNP_Scenario, Country_Depth_cm_CNP_Scenario) %>% # avoids lots of phyloseq::merge_samples warnings
  merge_samples2(group = "Country_Depth_cm_CNP_Scenario", funs = list()) %>%
  tax_glom(taxrank = "Phylum")

taxa_16s <- all_16s_phyloseq_rare_merged_phyla %>%
  tax_filter(
    tax_level = "Phylum", min_prevalence = 0.001,
    prev_detection_threshold = 100
  ) %>%
  comp_barplot("Phylum", n_taxa = 10, merge_other = TRUE, bar_width = 0.9, bar_outline_colour = NA, label = "Depth_cm", x = "Depth_cm",
               palette = taxa_pal_10) +
  facet_grid(CNP_Scenario ~ Country, scales = "free_x") +
  my_theme + theme(legend.text = element_text(face = "italic"), legend.position = "right") +
  ggtitle("A") +
  labs(x = "Depth (cm)", y = "Relative abundance", fill = "Prokaryotic Phylum")

# replot with correct taxa order
taxa_16s$layers[[1]] <- NULL
taxa_16s <- taxa_16s + geom_col(position = ggplot2::position_stack(reverse = TRUE))
taxa_16s
# manipulate data frame and save
taxa_16s$data %>%
  select(OTU, Abundance, Country, Depth_cm, CNP_Scenario) %>%
  dplyr::rename("Phylum" = OTU) %>%
  arrange(CNP_Scenario, Country, Depth_cm, desc(Abundance)) %>%
  writexl::write_xlsx(path = "/~your_name.xlsx")

# repeat for ITS
# "Phylum" can be also replaced with "Class", "Family", etc

