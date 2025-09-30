# effect of P rate on the frequency of occurrence of fungal taxa
# 0-10 cm
# transform
all_ITS_phyloseq_rare_df <- metagMisc::phyloseq_otu_occurrence(all_ITS_phyloseq_rare_FOO %>% ps_filter(Country == "Ireland", Depth_cm == "0-10"), variable = "CNP_Scenario", taxa_frequency = FALSE) %>%
  ps_melt() # adjust here
# KW
kw_ITS_0_10 <- all_ITS_phyloseq_rare_df %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  filter(n_distinct(Sample) > 1) %>%
  kruskal_test(Abundance ~ Sample) %>%
  adjust_pvalue(method = "fdr") %>% 
  filter(p.adj < 0.05)
kw_ITS_0_10
# filter for significantle affected genera
only_sig_ITS <- all_ITS_phyloseq_rare_df %>%
  filter(Genus %in% kw_ITS_0_10$Genus) %>%
  mutate(Sample = factor(Sample, levels = c("Low", "Medium", "High")))

# plot
sig_ITS_0_10 <- ggplot(only_sig_ITS, aes(y = Abundance, x = Sample)) +
  geom_bar(stat = "identity", aes(fill = Sample)) +
  scale_fill_manual(name = "", values = c(col2, "#D3AF37", col1)) +
  my_theme + theme(legend.position = "none") +
  labs(y = "FOO (%)", x = "P fertilization rate") +
  facet_wrap(~ Genus, scales = "free_y") +
  theme(strip.text.x = element_text(face = "italic", size = 8)) +
  stat_compare_means(method = "wilcox", comparisons = list(c("Low", "Medium"), c("Low", "High")), label = "p.signif", label.y = c(3.3, 3.6)) +
  ylim(0, 4)
sig_ITS_0_10
# manipulate data frames and save
sig_ITS_0_10$data %>%
  filter(Abundance > 0) %>%
  group_by(Sample) %>%
  summarize_at(vars(OTU), funs(n_distinct(.)))

sig_ITS_0_10$data %>%
  filter(Abundance > 0) %>%
  select(Sample, G_species, OTU, Abundance) %>%
  group_by(Sample, G_species) %>%
  filter(OTU == unique(OTU)) %>%
  mutate(Abundance = round(Abundance, 2)) %>%
  dplyr::rename("SP dose" = Sample, "Species" = G_species, "FOO (%)" = Abundance) %>%
  select(-OTU) %>%
  writexl::write_xlsx(path = "~/your_name.xlsx")