# alpha diversity plots comparing land uses
# these data were tested with shapiro wilcox test beforehand

depth_labeller <- as_labeller(c("0-10" = "0-10 cm", "10-30" = "10-30 cm", "30-50" = "30-50 cm"))
# 16S
# observed
shapiro.test(alpha_meta_iconica_16s$Observed)

significance <- alpha_meta_iconica_16s %>%
  group_by(Depth_cm) %>%
  t_test(Observed ~ Site_Type) %>%
  adjust_pvalue()
significance

stat_text <- data.frame(
  Depth_cm = c("0-10", "10-30", "30-50"),
  x = c(2.5, 2.5, 2),
  y = c(1400, 1400, 1400),
  label = c("p < 0.001", "ns", "ns")
)
alpha_site_16S <- ggplot(alpha_meta_iconica_16s, aes(x = Country_Code, y = Observed)) +
  geom_boxplot(aes(fill = Site_Type), outlier.shape = NA) +  
  ylab("Species richness") +
  xlab("Country") +
  my_theme + theme(legend.position = "none", axis.title.x = element_blank(), strip.text.x = element_text(size = 12)) +
  scale_color_manual(name = "Land type", labels = c("Arable", "Grassland"), values = c(col6, col7)) +
  scale_fill_manual(name = "Land type", labels = c("Arable", "Grassland"), values = c(col6, col7)) +
  geom_text(data = stat_text, aes(x = x, y = y, label = label), inherit.aes = FALSE, fontface = "italic", size = 3.5) +
  ggtitle("A") +
  facet_wrap(~ Depth_cm, labeller = depth_labeller, scales = "free_x") +
  ylim(200, 1400)

# shannon
shapiro.test(alpha_meta_iconica_16s$Shannon)

significance <- alpha_meta_iconica_16s %>%
  group_by(Depth_cm) %>%
  t_test(Shannon ~ Site_Type) %>%
  adjust_pvalue()
significance

stat_text <- data.frame(
  Depth_cm = c("0-10", "10-30", "30-50"),
  x = c(2.5, 2.5, 2),
  y = c(7, 7, 7),
  label = c("p < 0.001", "p < 0.001", "p < 0.001")
)
alpha_site_16S_shannon <- ggplot(alpha_meta_iconica_16s, aes(x = Country_Code, y = Shannon)) +
  geom_boxplot(aes(fill = Site_Type), outlier.shape = NA) +  
  ylab("Shannon diversity") +
  xlab("Country") +
  my_theme + theme(legend.position = "none", axis.title.x = element_blank(), strip.text.x = element_text(size = 12)) +
  scale_color_manual(name = "Land type", labels = c("Arable", "Grassland"), values = c(col6, col7)) +
  scale_fill_manual(name = "Land type", labels = c("Arable", "Grassland"), values = c(col6, col7)) +
  geom_text(data = stat_text, aes(x = x, y = y, label = label), inherit.aes = FALSE, fontface = "italic", size = 3.5) +
  ggtitle("B") +
  facet_wrap(~ Depth_cm, labeller = depth_labeller, scales = "free_x") +
  ylim(4.5, 7)


# ITS
# observed
shapiro.test(alpha_meta_iconica_ITS$Observed)

significance <- alpha_meta_iconica_ITS %>%
  group_by(Depth_cm) %>%
  kruskal_test(Observed ~ Site_Type) %>%
  adjust_pvalue()
significance

stat_text <- data.frame(
  Depth_cm = c("0-10", "10-30", "30-50"),
  x = c(2.5, 2.5, 2.5),
  y = c(300, 300, 300),
  label = c("p < 0.05", "p < 0.001", "p < 0.05")
)
alpha_site_ITS <- ggplot(alpha_meta_iconica_ITS, aes(x = Country_Code, y = Observed)) +
  geom_boxplot(aes(fill = Site_Type), outlier.shape = NA) +  
  ylab("Species richness") +
  xlab("Country") +
  my_theme + theme(legend.position = "none", axis.title.x = element_blank(), strip.text.x = element_text(size = 12)) +
  scale_color_manual(name = "Land type", labels = c("Arable", "Grassland"), values = c(col6, col7)) +
  scale_fill_manual(name = "Land type", labels = c("Arable", "Grassland"), values = c(col6, col7)) +
  geom_text(data = stat_text, aes(x = x, y = y, label = label), inherit.aes = FALSE, fontface = "italic", size = 3.5) +
  ggtitle("C") +
  ylim(0, 300) +
  facet_wrap(~ Depth_cm, labeller = depth_labeller, scales = "free_x")

# shannon
shapiro.test(alpha_meta_iconica_ITS$Shannon)

significance <- alpha_meta_iconica_ITS %>%
  group_by(Depth_cm) %>%
  kruskal_test(Shannon ~ Site_Type) %>%
  adjust_pvalue()
significance

stat_text <- data.frame(
  Depth_cm = c("0-10", "10-30", "30-50"),
  x = c(2.5, 2.5, 2.5),
  y = c(5, 5, 5),
  label = c("p < 0.001", "p < 0.001", "p < 0.001")
)
alpha_site_ITS_shannon <- ggplot(alpha_meta_iconica_ITS, aes(x = Country_Code, y = Shannon)) +
  geom_boxplot(aes(fill = Site_Type), outlier.shape = NA) +  
  ylab("Shannon diversity") +
  xlab("Country") +
  my_theme + theme(legend.position = "right", axis.title.x = element_blank(), strip.text.x = element_text(size = 12)) +
  scale_color_manual(name = "Land type", labels = c("Arable", "Grassland"), values = c(col6, col7)) +
  scale_fill_manual(name = "Land type", labels = c("Arable", "Grassland"), values = c(col6, col7)) +
  geom_text(data = stat_text, aes(x = x, y = y, label = label), inherit.aes = FALSE, fontface = "italic", size = 4) +
  ggtitle("D") +
  facet_wrap(~ Depth_cm, labeller = depth_labeller, scales = "free_x") +
  ylim(0, 5)
alpha_site_ITS_shannon