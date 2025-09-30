# alpha diversity
# 16S
alpha_summary_iconica_16s <- estimate_richness(all_16s_phyloseq_rare, measures = c("Observed", "Shannon"))

# this calculates eveness and adds it to the prior file

evenness_16s <- microbiome::evenness(all_16s_phyloseq_rare, "pielou")
alpha_summary_iconica_16s$Pielou <- evenness_16s$pielou

# combine results with sample_data. This is useful if we wish to plot later against other variables or information.
alpha_meta_iconica_16s <- data.frame(alpha_summary_iconica_16s, sample_data(all_16s_phyloseq_rare))

# convert data from wide to long for easier plotting
alpha_meta_iconica_16s_long = gather(alpha_meta_iconica_16s, Alpha, alpha_value, c("Observed", "Shannon", "Pielou"), factor_key = FALSE)

# summarise data
iconica_16s_observed_sum = summarySE(alpha_meta_iconica_16s_long, measurevar = "alpha_value", groupvars = c("Alpha", "Country", "Depth_cm", "CNP_Scenario", "Site_Type"))

iconica_16s_observed_sum$CNP_Scenario <- factor(iconica_16s_observed_sum$CNP_Scenario, levels = c("Low", "Medium", "High"))
# iconica_16s_observed_sum$P_dose <- factor(iconica_16s_observed_sum$P_dose, levels = c("0", "15.6", "23", "30", "41", "45", "105"))

alpha_meta_iconica_16s <- alpha_meta_iconica_16s %>%
  mutate(Country_Code = case_when(Country == "Denmark" ~ "DK", Country == "Sweden" ~ "SW", Country == "Ireland" ~ "IR", Country == "Netherlands" ~ "NL"),
         Country_Code = factor(Country_Code, levels = c("DK", "SW", "IR", "NL")))

# ITS
alpha_summary_iconica_ITS <- estimate_richness(all_ITS_phyloseq_rare, measures = c("Observed", "Shannon"))

# this calculates eveness and adds it to the prior file

evenness_ITS <- microbiome::evenness(all_ITS_phyloseq_rare, "pielou")
alpha_summary_iconica_ITS$Pielou <- evenness_ITS$pielou

# combine results with sample_data. This is useful if we wish to plot later against other variables or information.
alpha_meta_iconica_ITS <- data.frame(alpha_summary_iconica_ITS, sample_data(all_ITS_phyloseq_rare))

# convert data from wide to long for easier plotting
alpha_meta_iconica_ITS_long = gather(alpha_meta_iconica_ITS, Alpha, alpha_value, c("Observed", "Shannon", "Pielou"), factor_key = FALSE)

# summarise data
iconica_ITS_observed_sum = summarySE(alpha_meta_iconica_ITS_long, measurevar = "alpha_value", groupvars = c("Alpha", "Country", "Site_Type", "Depth_cm", "CNP_Scenario"))

iconica_ITS_observed_sum$CNP_Scenario <- factor(iconica_ITS_observed_sum$CNP_Scenario, levels = c("Low", "Medium", "High"))
# iconica_ITS_observed_sum$P_dose <- factor(iconica_ITS_observed_sum$P_dose, levels = c("0", "15.6", "23", "30", "41", "45", "105"))

alpha_meta_iconica_ITS <- alpha_meta_iconica_ITS %>%
  mutate(Country_Code = case_when(Country == "Denmark" ~ "DK", Country == "Sweden" ~ "SW", Country == "Ireland" ~ "IR", Country == "Netherlands" ~ "NL"),
         Country_Code = factor(Country_Code, levels = c("DK", "SW", "IR", "NL")))