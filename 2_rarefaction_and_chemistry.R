# rerefaction and cleaning
all_16s_phyloseq_rare <- rarefy_even_depth(all_16s_phyloseq, 12000, rngseed = TRUE)
all_16s_phyloseq_rare <- name_na_taxa(all_16s_phyloseq_rare, na_label = "Unknown <tax> (<rank>)")

all_16s_phyloseq_rare@sam_data$Site_Code = factor(all_16s_phyloseq_rare@sam_data$Site_Code, levels = c("JY", "LS1", "LS2", "JC1", "JC2", "LE"))
sample_data(all_16s_phyloseq_rare)$CNP_Scenario <- factor(sample_data(all_16s_phyloseq_rare)$CNP_Scenario, levels = c("Low", "Medium", "High"))


sample_data(all_16s_phyloseq_rare)$Country <- gsub(sample_data(all_16s_phyloseq_rare)$Country, pattern = "Denmark ", replacement = "Denmark")
sample_data(all_16s_phyloseq_rare)$Country <- gsub(sample_data(all_16s_phyloseq_rare)$Country, pattern = "Sweden ", replacement = "Sweden")
sample_data(all_16s_phyloseq_rare)$Country <- gsub(sample_data(all_16s_phyloseq_rare)$Country, pattern = "Ireland ", replacement = "Ireland")
sample_data(all_16s_phyloseq_rare)$Country <- gsub(sample_data(all_16s_phyloseq_rare)$Country, pattern = "Netherlands ", replacement = "Netherlands")
sample_data(all_16s_phyloseq_rare)$Country <- factor(sample_data(all_16s_phyloseq_rare)$Country, levels = c("Denmark", "Sweden", "Ireland", "Netherlands"))

all_ITS_phyloseq_rare <- rarefy_even_depth(all_ITS_phyloseq, 5000, rngseed = TRUE)
all_ITS_phyloseq_rare <- name_na_taxa(all_ITS_phyloseq_rare, na_label = "Unknown <tax> (<rank>)")
all_ITS_phyloseq_rare@sam_data$Site_Code = factor(all_ITS_phyloseq_rare@sam_data$Site_Code, levels = c("JY", "LS1", "LS2", "JC1", "JC2", "LE"))
sample_data(all_ITS_phyloseq_rare)$CNP_Scenario <- factor(sample_data(all_ITS_phyloseq_rare)$CNP_Scenario, levels = c("Low", "Medium", "High"))

all_16s_phyloseq_rare@sam_data$Ca_Mehlich3_mg.kg <- as.numeric(all_16s_phyloseq_rare@sam_data$Ca_Mehlich3_mg.kg)
all_16s_phyloseq_rare@sam_data$Soil_Organic_C_. <- as.numeric(all_16s_phyloseq_rare@sam_data$Soil_Organic_C_.)

all_ITS_phyloseq_rare@sam_data$Ca_Mehlich3_mg.kg <- as.numeric(all_ITS_phyloseq_rare@sam_data$Ca_Mehlich3_mg.kg)
all_ITS_phyloseq_rare@sam_data$Soil_Organic_C_. <- as.numeric(all_ITS_phyloseq_rare@sam_data$Soil_Organic_C_.)

sample_data(all_ITS_phyloseq_rare)$Country <- gsub(sample_data(all_ITS_phyloseq_rare)$Country, pattern = "Denmark ", replacement = "Denmark")
sample_data(all_ITS_phyloseq_rare)$Country <- gsub(sample_data(all_ITS_phyloseq_rare)$Country, pattern = "Sweden ", replacement = "Sweden")
sample_data(all_ITS_phyloseq_rare)$Country <- gsub(sample_data(all_ITS_phyloseq_rare)$Country, pattern = "Ireland ", replacement = "Ireland")
sample_data(all_ITS_phyloseq_rare)$Country <- gsub(sample_data(all_ITS_phyloseq_rare)$Country, pattern = "Netherlands ", replacement = "Netherlands")
sample_data(all_ITS_phyloseq_rare)$Country <- factor(sample_data(all_ITS_phyloseq_rare)$Country, levels = c("Denmark", "Sweden", "Ireland", "Netherlands"))

# creating metadata table
iconica_all_field_metadata <- data.frame(sample_data(all_ITS_phyloseq_rare))

# chemistry analysis
# SOC
metadata_sum_C <- iconica_all_field_metadata %>%
  group_by(Country, CNP_Scenario, Depth_cm) %>%
  dplyr::summarise(
    N = dplyr::n(),
    median = median(Soil_Organic_C_., na.rm = TRUE),
    IQR = IQR(Soil_Organic_C_., na.rm = TRUE),
    Nutrient = "SOC (%)",
    .groups = "drop"
  )

# TN
metadata_sum_N <- iconica_all_field_metadata %>%
  group_by(Country, CNP_Scenario, Depth_cm) %>%
  dplyr::summarise(
    N = dplyr::n(),
    median = median(Total_N_., na.rm = TRUE),
    IQR = IQR(Total_N_., na.rm = TRUE),
    Nutrient = "TN (%)",
    .groups = "drop"
  )

# mehlich 3 P
metadata_sum_P <- iconica_all_field_metadata %>%
  group_by(Country, CNP_Scenario, Depth_cm) %>%
  dplyr::summarise(
    N = dplyr::n(),
    median = median(P_Mehlich3_mg.kg, na.rm = TRUE),
    IQR = IQR(P_Mehlich3_mg.kg, na.rm = TRUE),
    Nutrient = "P (mg/kg)",
    .groups = "drop"
  )

metadata_sum <- bind_rows(metadata_sum_C, metadata_sum_N, metadata_sum_P)
# plot CNP
sum_plot <- ggplot(metadata_sum %>%
                     mutate(Depth_cm = factor(Depth_cm, levels = c("30-50", "10-30", "0-10")),
                            Nutrient = factor(Nutrient, levels = c("SOC (%)", "TN (%)", "P (mg/kg)"))),
                   aes(x = Depth_cm, y = median, colour = Country, group = Country, shape = Country)) +
  geom_point(size = 3) +
  geom_line(data = metadata_sum,
            aes(x = Depth_cm, y = median, group = Country, color = Country), linetype = "dashed") +
  geom_errorbar(aes(ymin = median - metadata_sum$IQR, 
                    ymax = median + metadata_sum$IQR), 
                width = 0.1) +
  my_theme +
  labs(x = "Depth (cm)", color = "Country", y = "") +
  facet_grid2(CNP_Scenario ~ Nutrient, scales = "free") +
  theme(legend.position = "right",
        axis.ticks = element_blank()) +
  coord_flip() +
  scale_color_manual(name = "Country", values = c(col2, col6, col3, col1)) +
  scale_shape_manual(name = "Country", values = c(17, 18, 15, 16))
sum_plot

# PCA
pca <- iconica_all_field_metadata %>%
  select(Sample_Code, Site_Type, Country, CNP_Scenario, Depth_cm, Sand_., Silt_., Clay_., pH, Soil_Organic_C_.,
         Total_N_., Al_Mehlich3_mg.kg, Ca_Mehlich3_mg.kg, Fe_Mehlich3_mg.kg, K_Mehlich3_mg.kg, Mg_Mehlich3_mg.kg, P_Mehlich3_mg.kg) %>%
  dplyr::rename("Clay" = Clay_., "Silt" = Silt_., "Sand" = Sand_., "C" = Soil_Organic_C_., "N" = Total_N_., "P" = P_Mehlich3_mg.kg,
                "Al" = Al_Mehlich3_mg.kg, "Ca" = Ca_Mehlich3_mg.kg, "Fe" = Fe_Mehlich3_mg.kg,
                "K" = K_Mehlich3_mg.kg, "Mg" = Mg_Mehlich3_mg.kg)


pca_codes <- prcomp(na.omit(pca[,c(6:17)]), scale = TRUE)
summary(pca_codes)
var_explained <- pca_codes$sdev^2/sum(pca_codes$sdev^2)
var_explained[1:5]

pca_codes_correlation <- data.frame(pca_codes$rotation[,1:2]) %>%
  rownames_to_column() %>%
  dplyr::rename("Variable" = rowname) %>%
  arrange(desc(PC1))
pca_codes_correlation

plot_pca <- autoplot(pca_codes,
                     data=pca,
                     colour="CNP_Scenario",
                     shape = "Depth_cm",
                     loadings = TRUE,
                     loadings.colour = 'black',
                     loadings.label.colour='black',
                     loadings.label = TRUE,
                     loadings.label.size = 4,
                     loadings.label.repel=TRUE,
                     size = 3) +
  my_theme + theme(legend.position = "right") +
  labs(color = "P fertilization\nrate") +
  scale_color_manual(name = "P fertilization\nrate", values = c(col2, "#D3AF37", col1), 
                     labels = c("Low", "Medium", "High")) +
  scale_shape_manual(name = "Depth (cm)", values = c(16, 17, 15), labels = c("0-10", "10-30", "30-50")) + 
  geom_mark_ellipse(aes(group = Country, linetype = Country, label = Country), alpha = 1, con.cap = 0.01, label.fontsize = 8,
                    show.legend = FALSE, con.type = "straight") +
  xlim(-0.25, 0.25) +
  ylim(-0.25, 0.25) 
plot_pca