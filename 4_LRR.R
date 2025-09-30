# combine 16S and ITS alpha diversity data
lrr_data <- alpha_meta_iconica_16s %>%
  select(Observed, Shannon, Sample_Code) %>%
  left_join(alpha_meta_iconica_ITS, by = "Sample_Code") %>%
  select(Sample_Code, Site_Type, Country_Code, CNP_Scenario, Depth_cm, Clay_., pH, Soil_Organic_C_., Total_N_., Al_Mehlich3_mg.kg, Ca_Mehlich3_mg.kg, Fe_Mehlich3_mg.kg, K_Mehlich3_mg.kg, Mg_Mehlich3_mg.kg, P_Mehlich3_mg.kg, Observed.x, Observed.y, Shannon.x, Shannon.y,) %>%
  dplyr::rename("Observed (16S)" = Observed.x, "Shannon (16S)" = Shannon.x,
                "Observed (ITS)" = Observed.y, "Shannon (ITS)" = Shannon.y,
                "Clay" = Clay_., "C" = Soil_Organic_C_., "N" = Total_N_., "P" = P_Mehlich3_mg.kg,
                "Al" = Al_Mehlich3_mg.kg, "Ca" = Ca_Mehlich3_mg.kg, "Fe" = Fe_Mehlich3_mg.kg,
                "K" = K_Mehlich3_mg.kg, "Mg" = Mg_Mehlich3_mg.kg) %>%
  filter(!is.na(pH)) %>%
  relocate(P, .before = Al)

data_gathered <- lrr_data %>%
  gather(Parameter, Measurement, Clay:"Shannon (ITS)", factor_key = TRUE)

# write function to calculate LRR, SE, and p-value
LRRd <- function(control, treatment) {
  # remove NAs
  control <- control[!is.na(control)]
  treatment <- treatment[!is.na(treatment)]
  
  # check if there’s enough data
  if (length(control) == 0 || length(treatment) == 0) 
    return(list(LRR = NA, SE = NA, p_value = NA))
  
  # check for non-positive values
  if (any(control <= 0) || any(treatment <= 0)) 
    return(list(LRR = NA, SE = NA, p_value = NA))
  
  # calculate means
  mean_control <- mean(control)
  mean_treatment <- mean(treatment)
  
  # calculate LRR
  LRR <- log(mean_treatment / mean_control)
  
  # calculate SE using the delta method
  sd_control <- sd(control)
  sd_treatment <- sd(treatment)
  n_control <- length(control)
  n_treatment <- length(treatment)
  
  SE <- sqrt((sd_treatment^2) / (n_treatment * mean_treatment^2) +
               (sd_control^2) / (n_control * mean_control^2))
  
  # calculate z-score and two-tailed p-value
  if (!is.na(SE) && SE != 0) {
    z <- LRR / SE
    p_value <- 2 * (1 - pnorm(abs(z)))
  } else {
    p_value <- NA
  }
  
  return(list(LRR = LRR, SE = SE, p_value = p_value))
}

LRR_Iconica <- function(df, country, treatment, depth, value) {
  Controls <- df %>%
    filter(Country_Code == country, CNP_Scenario == "Low", Depth_cm == depth, Parameter == value) %>%
    select(Measurement)
  
  Response <- df %>%
    filter(Country_Code == country, CNP_Scenario == treatment, Depth_cm == depth, Parameter == value) %>%
    select(Measurement)
  
  res <- LRRd(Controls$Measurement, Response$Measurement)
  return(res)
}

# loop to calculate LRRs
results_list <- list()

countries <- c("DK", "IR", "SW", "NL")
treatments <- c("Medium", "High")
depths <- c("0-10", "10-30", "30-50")
parameters <- c("Observed (16S)", "Shannon (16S)", "Observed (ITS)", "Shannon (ITS)")

row_index <- 1
for (country in countries) {
  for (treatment in treatments) {
    for (depth in depths) {
      for (parameter in parameters) {
        result <- tryCatch({
          LRR_value <- LRR_Iconica(df = data_gathered, 
                                   country = country, 
                                   treatment = treatment, 
                                   depth = depth, 
                                   value = parameter)
          data.frame(
            Country = country,
            Treatment = treatment,
            Depth = depth,
            Parameter = parameter,
            LRR = LRR_value$LRR,
            SE = LRR_value$SE,
            p_value = LRR_value$p_value
          )
        }, error = function(e) {
          data.frame(
            Country = country,
            Treatment = treatment,
            Depth = depth,
            Parameter = parameter,
            LRR = NA,
            SE = NA,
            p_value = NA
          )
        })
        results_list[[row_index]] <- result
        row_index <- row_index + 1
      }
    }
  }
}

# combine all results into a single dataframe
LRR_results_df <- do.call(rbind, results_list)
labeller <- as_labeller(c("Denmark" = "DK", "Sweden" = "SW", "Ireland" = "IR", "Netherlands" = "NL", "0-10" = "0-10 cm", "10-30" = "10-30 cm", "30-50" = "30-50 cm", "Arable" = "Arable", "Grassland" = "Grassland"))

LRR_results_df <- LRR_results_df %>%
  mutate(Site_Type = case_when(Country %in% c("DK", "SW") ~ "Arable", TRUE ~ "Grassland"),
         Parameter = factor(Parameter, levels = c("Shannon (ITS)", "Observed (ITS)", "Shannon (16S)", "Observed (16S)")),
         Site_Type = factor(Site_Type, levels = c("Grassland", "Arable")),
         Country = factor(Country, levels = c("DK", "SW", "IR", "NL")))
# plot LRR
LRR_plot_SP <- ggplot(LRR_results_df, aes(Parameter, LRR, colour = Treatment)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_point(size = 3, aes(shape = Country)) + geom_linerange(aes(ymin = LRR - SE, ymax = LRR + SE)) + 
  coord_flip() + ylab("Log R ratio") + xlab("Alpha diversity") +
  my_theme +
  scale_shape_manual(name = "Country", values = c(17, 18, 15, 16)) +
  scale_color_manual(name = "P fertilization\nrate", values = c("#D3AF37", col1), labels = c("Medium", "High")) +
  facet_grid(Depth ~ Site_Type, scales = "free", labeller = labeller)

LRR_plot_SP