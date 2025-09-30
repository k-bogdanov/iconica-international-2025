# set seed
set.seed(1509)

# ordinate 16S
all_16s_NMDS_ord <- ordinate(all_16s_phyloseq_rare, "NMDS", "bray")
all_16s_NMDS_ord

stressplot(all_16s_NMDS_ord)

# P fertilization rate
all_16s_phyloseq_rare@sam_data$desc_P <- c("Depth: p ≤ 0.001, R² = 0.406ns\nP rate: ns",
                                           "Depth: ns\nP rate: ns",
                                           "Depth: p ≤ 0.001, R² = 0.289\nP rate: ns",
                                           "Depth: ns\nP rate: ns")[all_16s_phyloseq_rare@sam_data$Country] # from adonis tests
# create global label
nmds_label_16s = grobTree(textGrob("Global stats:\nCountry: p ≤ 0.001, R² = 0.391\nLand use: p ≤ 0.001, R² = 0.175\nDepth: p ≤ 0.001, R² = 0.089\nP rate: ns", x = 0.01,  y = 0.90, hjust = 0, gp = gpar(col = "Black", fontsize = 8, fontface = "plain")))

nmds_plot_16s <- plot_ordination(all_16s_phyloseq_rare, all_16s_NMDS_ord, justDF = TRUE) %>% 
  ggplot(aes(NMDS1, NMDS2)) +
  scale_color_manual(name = "P fertilization\nrate", values = c(col2, "#D3AF37", col1), labels = c("Low", "Medium", "High")) +
  scale_shape_manual(name = "Depth (cm)", values = c(16, 17, 15), labels = c("0-10", "10-30", "30-50")) +  
  geom_point(aes(color = CNP_Scenario, shape = Depth_cm), size = 3) +
  my_theme + theme(legend.position = "none") +
  annotation_custom(nmds_label_16s) +
  geom_mark_ellipse(aes(group = Country, linetype = Country, label = Country, description = desc_P), alpha = 1, con.cap = 0.01, label.fontsize = 8, show.legend = FALSE, con.type = "straight") +
  ggtitle("A") +
  xlim(-5.5, 5.5) +
  ylim(-5.5, 5.5)
nmds_plot_16s

# SOC %
all_16s_phyloseq_rare@sam_data$desc_SOC <- c("SOC: p = 0.002, R² = 0.174\nDepth: p ≤ 0.001, R² = 0.406",
                                             "SOC: p = 0.006, R² = 0.121\nDepth: ns",
                                             "SOC: p ≤ 0.001, R² = 0.201\nDepth: p ≤ 0.001, R² = 0.289",
                                             "SOC: p = 0.015, R² = 0.337\nDepth: ns")[all_16s_phyloseq_rare@sam_data$Country] 

nmds_plot_16s_SOC <- plot_ordination(all_16s_phyloseq_rare, all_16s_NMDS_ord, justDF = TRUE) %>% 
  ggplot(aes(NMDS1, NMDS2)) +
  scale_color_wa_c(name = "SOC (%)", reverse = F, palette = "vantage") +
  scale_shape_manual(name = "Depth (cm)", values = c(16, 17, 15), labels = c("0-10", "10-30", "30-50")) +  
  geom_point(aes(color = Soil_Organic_C_., shape = Depth_cm), alpha = 1, size = 3) +
  my_theme +
  geom_mark_ellipse(aes(group = Country, linetype = Country, label = Country, description = desc_SOC), alpha = 1, con.cap = 0.01, label.fontsize = 8,   show.legend = FALSE, con.type = "straight") +
  labs(color = "Land type") +
  ggtitle("Prokaryotes") +
  xlim(-5.5, 5.5) +
  ylim(-5.5, 5.5)
nmds_plot_16s_SOC

# TN %
all_16s_phyloseq_rare@sam_data$desc_TN<- c("TN: p = 0.002, R² = 0.167\nDepth: p ≤ 0.001, R² = 0.406",
                                           "TN: p = 0.044, R² = 0.104\nDepth: ns",
                                           "TN: p ≤ 0.001, R² = 0.201\nDepth: p ≤ 0.001, R² = 0.289",
                                           "TN: p = 0.020, R² = 0.356\nDepth: ns")[all_16s_phyloseq_rare@sam_data$Country]

nmds_plot_16s_TN <- plot_ordination(all_16s_phyloseq_rare, all_16s_NMDS_ord, justDF = TRUE) %>% 
  ggplot(aes(NMDS1, NMDS2)) +
  scale_color_wa_c(name = "TN (%)", reverse = T, palette = "ferries") +
  scale_shape_manual(name = "Depth (cm)", values = c(16, 17, 15), labels = c("0-10", "10-30", "30-50")) +  
  geom_point(aes(color = Total_N_., shape = Depth_cm), alpha = 1, size = 3) +
  my_theme +
  geom_mark_ellipse(aes(group = Country, linetype = Country, label = Country, description = desc_TN), alpha = 1, con.cap = 0.01, label.fontsize = 8,   show.legend = FALSE, con.type = "straight") +
  labs(color = "Land type") +
  ggtitle("Prokaryotes") +
  xlim(-5.5, 5.5) +
  ylim(-5.5, 5.5)
nmds_plot_16s_TN

# ordinate ITS
all_ITS_NMDS_ord <- ordinate(all_ITS_phyloseq_rare, "NMDS", "bray")
all_ITS_NMDS_ord

stressplot(all_ITS_NMDS_ord)

# P fertilization rate
all_ITS_phyloseq_rare@sam_data$desc_P <- c("Depth: p ≤ 0.001, R² = 0.222\nP rate: ns",
                                           "Depth: ns\nP rate: ns",
                                           "Depth: p ≤ 0.001, R² = 0.152\nP rate: p = 0.011, R² = 0.068",
                                           "Depth: ns\nP rate: ns")[all_ITS_phyloseq_rare@sam_data$Country]

# create global label
nmds_label_ITS = grobTree(textGrob("Global stats:\nCountry: p ≤ 0.001, R² = 0.372\nLand use: p ≤ 0.001, R² = 0.196\nDepth: p = 0.005, R² = 0.042\nP rate: p = 0.042, R² = 0.029", x = 0.01,  y = 0.1, hjust = 0, gp = gpar(col = "Black", fontsize = 8, fontface = "plain")))

nmds_plot_ITS <- plot_ordination(all_ITS_phyloseq_rare, all_ITS_NMDS_ord, justDF = TRUE) %>% 
  ggplot(aes(NMDS1, NMDS2)) +
  scale_color_manual(name = "P fertilization\nrate", values = c(col2, "#D3AF37", col1), labels = c("Low", "Medium", "High")) +
  scale_shape_manual(name = "Depth (cm)", values = c(16, 17, 15), labels = c("0-10", "10-30", "30-50")) +  
  geom_point(aes(color = CNP_Scenario, shape = Depth_cm), size = 3) +
  my_theme +
  annotation_custom(nmds_label_ITS) +
  geom_mark_ellipse(aes(group = Country, linetype = Country, label = Country, description = desc_P), alpha = 1, con.cap = 0.01, label.fontsize = 8, show.legend = FALSE, con.type = "straight") +
  ggtitle("B") +
  xlim(-5.5, 5.5) +
  ylim(-5.5, 5.5)
nmds_plot_ITS


# SOC %
all_ITS_phyloseq_rare@sam_data$desc_SOC <- c("SOC: ns\nDepth: p ≤ 0.001, R² = 0.222",
                                             "SOC: ns\nDepth: ns",
                                             "SOC: p ≤ 0.001, R² = 0.119\nDepth: p ≤ 0.001, R² = 0.152",
                                             "SOC: p = 0.013, R² = 0.296\nDepth: ns")[all_ITS_phyloseq_rare@sam_data$Country]

nmds_plot_ITS_SOC <- plot_ordination(all_ITS_phyloseq_rare, all_ITS_NMDS_ord, justDF = TRUE) %>% 
  ggplot(aes(NMDS1, NMDS2)) +
  scale_color_wa_c(name = "SOC (%)", reverse = F, palette = "vantage") +
  scale_shape_manual(name = "Depth (cm)", values = c(16, 17, 15), labels = c("0-10", "10-30", "30-50")) +  
  geom_point(aes(color = Soil_Organic_C_., shape = Depth_cm), alpha = 1, size = 3) +
  my_theme +
  geom_mark_ellipse(aes(group = Country, linetype = Country, label = Country, description = desc_SOC), alpha = 1, con.cap = 0.01, label.fontsize = 8,   show.legend = FALSE, con.type = "straight") +
  labs(color = "Land type") +
  ggtitle("Fungi") +
  xlim(-5.5, 5.5) +
  ylim(-5.5, 5.5)
nmds_plot_ITS_SOC



# TN %
all_ITS_phyloseq_rare@sam_data$desc_TN <- c("TN: p = 0.049, R² = 0.062\nDepth: p ≤ 0.001, R² = 0.222",
                                            "TN: ns\nDepth: ns",
                                            "TN: p ≤ 0.001, R² = 0.123\nDepth: p ≤ 0.001, R² = 0.152",
                                            "TN: p = 0.014, R² = 0.294\nDepth: ns")[all_ITS_phyloseq_rare@sam_data$Country]

nmds_plot_ITS_TN <- plot_ordination(all_ITS_phyloseq_rare, all_ITS_NMDS_ord, justDF = TRUE) %>% 
  ggplot(aes(NMDS1, NMDS2)) +
  scale_color_wa_c(name = "TN (%)", reverse = T, palette = "ferries") +
  scale_shape_manual(name = "Depth (cm)", values = c(16, 17, 15), labels = c("0-10", "10-30", "30-50")) +  
  geom_point(aes(color = Total_N_., shape = Depth_cm), alpha = 1, size = 3) +
  my_theme +
  geom_mark_ellipse(aes(group = Country, linetype = Country, label = Country, description = desc_TN), alpha = 1, con.cap = 0.01, label.fontsize = 8,   show.legend = FALSE, con.type = "straight") +
  labs(color = "Land type") +
  ggtitle("Fungi") +
  xlim(-5.5, 5.5) +
  ylim(-5.5, 5.5)
nmds_plot_ITS_TN

# then combine using patchwork
