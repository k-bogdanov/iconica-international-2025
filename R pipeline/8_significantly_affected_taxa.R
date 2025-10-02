# significantly affected 16S taxa
# only 16S, replace "16S" with "ITS"
# transform ps
all_16S_phyloseq_rare_df <- all_16s_phyloseq_rare %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 
# run KW test
sig_affected_phyla_16S <- all_16S_phyloseq_rare_df %>% 
  group_by(Country, Phylum) %>% 
  kruskal_test(Abundance ~ CNP_Scenario) %>% 
  adjust_pvalue(method = "fdr")%>% 
  filter(p.adj < 0.05)
sig_affected_phyla_16S
# filter for significantly affected and exclude zero-abundant
only_sig_16S <- all_16S_phyloseq_rare_df %>%
  filter(Phylum %in% sig_affected_phyla_16S$Phylum, Abundance > 0) # filtering at this stage ensures that we keep taxa which were not present
                                                                   # under any of P rates. this is important because we can see what
                                                                   # microorganisms were completely supressed/promoted by P

# calculate differences within treatments
low_abundance <- only_sig_16S %>%
  filter(CNP_Scenario == "Low") %>%
  select(OTU, Country, Depth_cm, Abundance, Phylum) %>%
  mutate(Abundance = Abundance * 100) %>%
  dplyr::rename("Abundance_Low" = Abundance)
low_abundance

medium_abundance <- only_sig_16S %>%
  filter(CNP_Scenario == "Medium") %>%
  select(OTU, Country, Depth_cm, Abundance, Phylum) %>%
  mutate(Abundance = Abundance * 100) %>%
  dplyr::rename("Abundance_Medium" = Abundance)
medium_abundance

high_abundance <- only_sig_16S %>%
  filter(CNP_Scenario == "High") %>%
  select(OTU, Country, Depth_cm, Abundance, Phylum) %>%
  mutate(Abundance = Abundance * 100) %>%
  dplyr::rename("Abundance_High" = Abundance)
high_abundance


# differences between treatments
diff_abundance_high <- low_abundance %>%
  left_join(high_abundance, by = c("OTU", "Country", "Depth_cm")) %>%
  mutate(Difference = Abundance_Low - Abundance_High, Diff_Group = "Abundance_High") %>%
  select(-Abundance_Low, -Abundance_High)

diff_abundance_medium <- low_abundance %>%
  left_join(medium_abundance, by = c("OTU", "Country", "Depth_cm")) %>%
  mutate(Difference = Abundance_Low - Abundance_Medium, Diff_Group = "Abundance_Medium") %>%
  select(-Abundance_Low, -Abundance_Medium)

diff_abundance_total <- bind_rows(diff_abundance_high, diff_abundance_medium)

abundances <- low_abundance %>%
  left_join(medium_abundance, by = c("OTU", "Country", "Depth_cm", "Phylum")) %>%
  left_join(high_abundance, by = c("OTU", "Country", "Depth_cm", "Phylum")) %>%
  relocate(Phylum, .before = Abundance_Low) %>%
  gather(CNP_Scenario, Abundance, Abundance_Low:Abundance_High, factor_key = TRUE)

# do Wilcoxon tests
phyla_in_all_groups <- abundances %>%
  filter(Country != "Sweden") %>%
  group_by(Country, Depth_cm, CNP_Scenario) %>%
  summarise(Phylum_list = list(unique(Phylum)), .groups = "drop") %>%
  summarise(common_phyla = reduce(Phylum_list, intersect)) %>%
  pull(common_phyla)

safe_wilcox_tests <- function(data) {
  comparisons <- list(
    c("Abundance_Low", "Abundance_High"),
    c("Abundance_Low", "Abundance_Medium")
  )
  
  map_dfr(comparisons, function(comp) {
    tryCatch({
      wilcox_test(data, Abundance ~ CNP_Scenario, 
                  comparisons = list(comp)) %>%
        mutate(group1 = comp[1], group2 = comp[2])
    }, error = function(e) {
      tibble()
    })
  })
}

results <- abundances %>%
  filter(Country != "Sweden", Phylum %in% phyla_in_all_groups, !is.na(Phylum)) %>%
  group_by(Phylum, Country, Depth_cm) %>%
  group_modify(~ {
    needed_levels <- c("Abundance_Low", "Abundance_Medium", "Abundance_High")
    if (all(c("Abundance_Low") %in% .x$CNP_Scenario) &&
        sum(.x$CNP_Scenario %in% needed_levels) >= 2) {
      .x <- .x %>% filter(CNP_Scenario %in% needed_levels)
      safe_wilcox_tests(.x)
    } else {
      tibble()
    }
  }) %>%
  ungroup() %>%
  adjust_pvalue(method = "fdr") %>%
  filter(p.adj < 0.05)

results

results_for_plot <- results %>%
  mutate(
    y.position = 1,  # adjust per group if needed
    p.label = ifelse(p.adj < 0.001, "***",
                     ifelse(p.adj < 0.01, "**",
                            ifelse(p.adj < 0.05, "*", "ns"))),
    Diff_Group = case_when(
      group1 == "Abundance_Low" & group2 == "Abundance_High" ~ "Abundance_High",
      group1 == "Abundance_Low" & group2 == "Abundance_Medium" ~ "Abundance_Medium",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::rename(
    Phylum.x = Phylum
  )

results_for_plot$Diff_Group <- factor(results_for_plot$Diff_Group, levels = c("Abundance_Medium", "Abundance_High"))
diff_abundance_total$Diff_Group <- factor(diff_abundance_total$Diff_Group, levels = c("Abundance_Medium", "Abundance_High"))

# plot countries individually
# denmark plot
diff_abundance_total_DK <- diff_abundance_total %>%
  filter(Country == "Denmark") %>%
  filter(Phylum.x != "NA" & Phylum.y != "NA" & Difference != "NA")

phylum_levels_DK <- levels(factor(diff_abundance_total_DK$Phylum.x))

results_for_plot_DK <- results_for_plot %>%
  filter(Country == "Denmark") %>%
  mutate(
    y.position = case_when(
      Depth_cm == "0-10" ~ 3,
      Depth_cm == "10-30" ~ 3,
      Depth_cm == "30-50" ~ 10
    ),
    Phylum.x = factor(Phylum.x, levels = phylum_levels_DK),
    x_numeric = as.numeric(Phylum.x),
    x_adj = case_when(
      Diff_Group == "Abundance_High" ~ x_numeric + 0.1,
      Diff_Group == "Abundance_Medium" ~ x_numeric - 0.3,
      TRUE ~ x_numeric
    )
  )
labeller <- as_labeller(c("Denmark" = "DK", "Sweden" = "SW", "Ireland" = "IR", "Netherlands" = "NL", "0-10" = "0-10 cm", "10-30" = "10-30 cm", "30-50" = "30-50 cm", "Arable" = "Arable", "Grassland" = "Grassland"))

difference_16S_DK <- ggplot(diff_abundance_total_DK, aes(x = Phylum.x, y = Difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_boxplot(aes(color = Diff_Group), outlier.alpha = 0.3) +
  facet_grid(~ Depth_cm, scales = "free_x", labeller = labeller) +
  geom_text(
    data = results_for_plot_DK,
    aes(x = x_adj, y = y.position, label = p.label, color = Diff_Group),
    inherit.aes = FALSE,
    size = 5,
    show.legend = FALSE) +
  my_theme + theme(axis.text.y = element_text(face = "italic"),
                   legend.position = "none") +
  scale_color_manual(name = "P fertilization\nrate", values = c("#D3AF37", col1), 
                     labels = c("Medium", "High")) +
  scale_fill_manual(name = "P fertilization\nrate", values = c("#D3AF37", col1), 
                    labels = c("Medium", "High")) +
  coord_flip() +
  labs(x = "", y = "") +
  ggtitle("DK (arable)")
difference_16S_DK

difference_16S_DK$data

# Ireland plot
diff_abundance_total_IR <- diff_abundance_total %>%
  filter(Country == "Ireland") %>%
  filter(Phylum.x != "NA" & Phylum.y != "NA" & Difference != "NA")

phylum_levels_IR <- levels(factor(diff_abundance_total_IR$Phylum.x))

results_for_plot_IR <- results_for_plot %>%
  filter(Country == "Ireland") %>%
  mutate(
    y.position = case_when(
      Depth_cm == "0-10" ~ 3,
      Depth_cm == "10-30" ~ 2,
      Depth_cm == "30-50" ~ 5
    ),
    Phylum.x = factor(Phylum.x, levels = phylum_levels_IR),
    x_numeric = as.numeric(Phylum.x),
    x_adj = case_when(
      Diff_Group == "Abundance_High" ~ x_numeric + 0.1,
      Diff_Group == "Abundance_Medium" ~ x_numeric - 0.3,
      TRUE ~ x_numeric
    )
  )

difference_16S_IR <- ggplot(diff_abundance_total_IR, aes(x = Phylum.x, y = Difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_boxplot(aes(color = Diff_Group), outlier.alpha = 0.3) +
  facet_grid(~ Depth_cm, labeller = labeller, scales = "free") +
  geom_text(
    data = results_for_plot_IR,
    aes(x = x_adj, y = y.position, label = p.label, color = Diff_Group),
    inherit.aes = FALSE,
    size = 5,
    show.legend = FALSE) +
  my_theme + theme(axis.text.y = element_text(face = "italic"),
                   legend.position = "none") +
  scale_color_manual(name = "P fertilization\nrate", values = c("#D3AF37", col1), 
                     labels = c("Medium", "High")) +
  scale_fill_manual(name = "P fertilization\nrate", values = c("#D3AF37", col1), 
                    labels = c("Medium", "High")) +
  coord_flip() +
  labs(x = "", y = "") +
  ggtitle("IR (grassland)")
difference_16S_IR

# Netherlands
diff_abundance_total_NL <- diff_abundance_total %>%
  filter(Country == "Netherlands") %>%
  filter(Phylum.x != "NA" & Phylum.y != "NA" & Difference != "NA")

phylum_levels_NL <- levels(factor(diff_abundance_total_NL$Phylum.x))

results_for_plot_NL <- results_for_plot %>%
  filter(Country == "Netherlands") %>%
  mutate(
    y.position = case_when(
      Depth_cm == "0-10" ~ 2,
      Depth_cm == "10-30" ~ 5,
      Depth_cm == "30-50" ~ -10
    ),
    Phylum.x = factor(Phylum.x, levels = phylum_levels_NL),
    x_numeric = as.numeric(Phylum.x),
    x_adj = case_when(
      Diff_Group == "Abundance_High" ~ x_numeric + 0.1,
      Diff_Group == "Abundance_Medium" ~ x_numeric - 0.3,
      TRUE ~ x_numeric
    )
  )

difference_16S_NL <- ggplot(diff_abundance_total_NL, aes(x = Phylum.x, y = Difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_boxplot(aes(color = Diff_Group), outlier.alpha = 0.3) +
  facet_grid(~ Depth_cm, labeller = labeller, scales = "free") +
  geom_text(
    data = results_for_plot_NL,
    aes(x = x_adj, y = y.position, label = p.label, color = Diff_Group),
    inherit.aes = FALSE,
    size = 5,
    show.legend = FALSE) +
  my_theme + theme(axis.text.y = element_text(face = "italic"),
                   legend.position = "right") +
  scale_color_manual(name = "P fertilization\nrate", values = c("#D3AF37", col1), 
                     labels = c("Medium", "High")) +
  scale_fill_manual(name = "P fertilization\nrate", values = c("#D3AF37", col1), 
                    labels = c("Medium", "High")) +
  coord_flip() +
  labs(x = "", y = "") +
  ggtitle("NL (grassland)")
difference_16S_NL

# Sweden
abundances_sweden <- abundances %>% filter(Country == "Sweden")
# do stats
phyla_in_all_groups <- abundances_sweden %>%
  group_by(Country, Depth_cm, CNP_Scenario) %>%
  summarise(Phylum_list = list(unique(Phylum)), .groups = "drop") %>%
  summarise(common_phyla = reduce(Phylum_list, intersect)) %>%
  pull(common_phyla)

safe_wilcox_tests <- function(data) {
  comparisons <- list(
    c("Abundance_Low", "Abundance_High"),
    c("Abundance_Low", "Abundance_Medium")
  )
  
  map_dfr(comparisons, function(comp) {
    tryCatch({
      wilcox_test(data, Abundance ~ CNP_Scenario, 
                  comparisons = list(comp)) %>%
        mutate(group1 = comp[1], group2 = comp[2])
    }, error = function(e) {
      tibble()
    })
  })
}

results_sweden <- abundances_sweden %>%
  filter(Phylum %in% phyla_in_all_groups, !is.na(Phylum)) %>%
  group_by(Phylum, Country, Depth_cm) %>%
  group_modify(~ {
    needed_levels <- c("Abundance_Low", "Abundance_Medium", "Abundance_High")
    if (all(c("Abundance_Low") %in% .x$CNP_Scenario) &&
        sum(.x$CNP_Scenario %in% needed_levels) >= 2) {
      .x <- .x %>% filter(CNP_Scenario %in% needed_levels)
      safe_wilcox_tests(.x)
    } else {
      tibble()
    }
  }) %>%
  ungroup() %>%
  adjust_pvalue(method = "fdr") %>%
  filter(p.adj < 0.05)


diff_abundance_total_SW <- diff_abundance_total %>%
  filter(Country == "Sweden") %>%
  filter(Phylum.x != "NA" & Phylum.y != "NA" & Difference != "NA")
phylum_levels_SW <- levels(factor(diff_abundance_total_SW$Phylum.x))

# make stats table for plotting
results_for_plot_SW <- results_sweden %>%
  mutate(
    Phylum = factor(Phylum, levels = phylum_levels_SW),
    x_numeric = as.numeric(Phylum),
    y.position = 2, 
    p.label = ifelse(p.adj < 0.001, "***",
                     ifelse(p.adj < 0.01, "**",
                            ifelse(p.adj < 0.05, "*", "ns"))),
    Diff_Group = case_when(
      group1 == "Abundance_Low" & group2 == "Abundance_High" ~ "Abundance_High",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::rename(
    Phylum.x = Phylum
  )

# plot
difference_16S_SW <- ggplot(diff_abundance_total_SW, aes(x = Phylum.x, y = Difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_boxplot(aes(color = Diff_Group), outlier.alpha = 0.3) +
  facet_grid(~ Depth_cm, labeller = labeller, scales = "free") +
  my_theme + theme(axis.text.y = element_text(face = "italic"),
                   legend.position = "none") +
  geom_text(
    data = results_for_plot_SW,
    aes(x = x_numeric, y = y.position, label = p.label, color = Diff_Group),
    inherit.aes = FALSE,
    size = 5,
    show.legend = FALSE) +
  scale_color_manual(name = "P fertilization\nrate", values = c(col1), 
                     labels = c("Medium", "High")) +
  scale_fill_manual(name = "P fertilization\nrate", values = c(col1), 
                    labels = c("Medium", "High")) +
  coord_flip() +
  labs(x = "", y = "") +
  ggtitle("SW (arable)")
difference_16S_SW
