# 1
# adonis (global)
# create a data frame using your sample_data
df_ado = as(sample_data(all_16s_phyloseq_rare %>% ps_filter(Country == "Ireland")), "data.frame") # adjust ps (16S or ITS) and other variables
# calculate your Bray distance matrix
site_ado = phyloseq::distance(all_16s_phyloseq_rare %>% ps_filter(Country == "Ireland"), "bray")
# perform your ADONIS test
site_ado_stats <- adonis2(site_ado ~ Depth_cm, df_ado)
# check results
site_ado_stats


# 2
# adonis (for each depth, loop)
all_16s_phyloseq_rare_ado <- all_16s_phyloseq_rare %>%
  ps_select(-P_Tot_Digest_mg.kg, -Al_Tot_Digest_mg.kg, -Ca_Tot_Digest_mg.kg, -Fe_Tot_Digest_mg.kg, -K_Tot_Digest_mg.kg, -Mg_Tot_Digest_mg.kg, -SOC_stocks_tC.ha)

# Ireland
# define depths
depth_levels <- c("0-10", "10-30", "30-50")

# start the loop
adonis_results <- list()

for (depth in depth_levels) {
  
  # filter and prepare data
  ps_subset <- all_16s_phyloseq_rare_ado %>%
    ps_filter(Country == "Ireland" & Depth_cm == depth)
  
  distance <- phyloseq::distance(ps_subset, "bray")
  
  ado_df <- as(sample_data(ps_subset %>%
                             ps_mutate(across(c(14:30), round, 1)) %>%
                             ps_select(c(10, 14:30))), "data.frame")
  
  factors <- colnames(ado_df)
  
  p_values <- numeric(length(factors))
  r_values <- numeric(length(factors))
  
  # loop over variables for adonis2
  for (i in seq_along(factors)) {
    ado_res <- adonis2(distance ~ ado_df[[factors[i]]], ado_df)
    p_values[i] <- ado_res$`Pr(>F)`[1]
    r_values[i] <- ado_res$R2[1]
  }
  
  # store results
  adonis_results[[depth]] <- data.frame(
    factors = factors,
    p_values = p_values,
    r_values = r_values,
    Country = "Ireland",
    Depth = depth
  )
}

# combine all depth-level results
Ireland_16S <- bind_rows(adonis_results) %>%
  mutate(Type = "16S")

Ireland_16S

# run for all countries and then combine with bind_rows()

# 3
# anosim
Country_group = get_variable(all_ITS_phyloseq_rare, "Depth_cm") # adjust here 16S or ITS and variable
Country_ano = anosim(phyloseq::distance(all_ITS_phyloseq_rare, "bray"), Country_group)
Country_ano$signif
Country_ano$statistic

?anosim