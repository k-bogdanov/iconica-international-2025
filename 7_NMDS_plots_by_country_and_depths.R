# NMDS plots split by country and depth
# 16S only, replace "16S" with "ITS" 
countries <- c("Denmark", "Sweden", "Ireland", "Netherlands")

depths_by_country <- list(
  Denmark = c("0-10", "10-30", "30-50"),
  Sweden = c("0-10", "10-30"),
  Ireland = c("0-10", "10-30", "30-50"),
  Netherlands = c("0-10", "10-30", "30-50")
)

nmds_plots <- list()

for (country in countries) {
  for (depth in depths_by_country[[country]]) {
    
    # define abbreviation for each country
    country_abbr <- dplyr::recode(country,
                                  "Ireland" = "IR (grassland)",
                                  "Denmark" = "DK (arable)",
                                  "Netherlands" = "NL (grassland)",
                                  "Sweden" = "SW (arable)")
    
    # filter phyloseq object
    ps_subset <- all_16s_phyloseq_rare %>% ps_filter(Country == country & Depth_cm == depth)
    
    # run NMDS ordination
    ord <- ordinate(ps_subset, "NMDS", "bray")
    
    # extract ordination data
    ord_df <- plot_ordination(ps_subset, ord, justDF = TRUE)
    
    # remove extreme NMDS1 point only for Denmark 30-50 cm
    if (country == "Denmark" && depth == "30-50") {
      # define a threshold to exclude extreme NMDS1
      ord_df <- ord_df %>% filter(NMDS1 < 5)
    }
    
    # get adonis stats label for this subset
    adonis_label <- adonis_16S %>% # it is important to have adonis done by this stage, otherwise won't work
      dplyr::rename("Depth_cm" = Depth) %>%
      filter(factors == "CNP_Scenario", Country == country, Depth_cm == depth, !is.na(p_values)) %>%
      mutate(across(c(p_values, r_values), round, 3),
             label = paste0("p = ", p_values, "\n", "R² = ", r_values)) %>%
      left_join(ord_df, by = c("Country", "Depth_cm")) %>%
      group_by(Country, Depth_cm, label) %>%
      summarise_at(vars(NMDS1, NMDS2), min)
    
    if (depth == "0-10") {
      plot_title <- country_abbr
    } else {
      plot_title <- NULL
    }
    
    # plot
    plot <- ggplot(ord_df, aes(NMDS1, NMDS2)) +
      geom_point(aes(color = CNP_Scenario), size = 3) +
      geom_text(data = adonis_label, aes(x = NMDS1, y = NMDS2, label = label),
                inherit.aes = FALSE, hjust = 0, vjust = 0, size = 3.2, fontface = "italic") +
      scale_color_manual(name = "P fertilization\nrate", values = c(col2, "#D3AF37", col1),
                         labels = c("Low", "Medium", "High")) +
      my_theme +
      theme(legend.position = "none") +
      labs(title = plot_title, x = "", y = "")
    
    # store in a list
    plot_key <- paste0(country, "_", gsub("-", "_", depth))
    nmds_plots[[plot_key]] <- plot
  }
}


plot_with_legend <- nmds_plots[[1]] + theme(legend.position = "right")
legend_plot <- get_legend(plot_with_legend)
legend_patch <- cowplot::plot_grid(legend_plot)

plot_list <- unname(nmds_plots)

plot_list_with_legend <- append(plot_list, list(legend_patch), after = 5)
# add depths
plot_list_with_legend[[7]] <- plot_list_with_legend[[7]] + ylab("0-10 cm")
plot_list_with_legend[[8]]  <- plot_list_with_legend[[8]] + ylab("10-30 cm")
plot_list_with_legend[[9]]  <- plot_list_with_legend[[9]] + ylab("30-50 cm")

design1 =
  "GJAD
 HKBE
 ILCF"
# combine with patchwork
combined_plot <- wrap_plots(plot_list_with_legend, ncol = 4, nrow = 3, design = design1) & theme(axis.text = element_blank())
combined_plot

