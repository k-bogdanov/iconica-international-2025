# theme
my_theme <- theme(axis.title.y = element_text(size = 12, face = "bold", vjust = 3),
                  axis.title.x = element_text(size = 12, face = "bold", vjust = -1),
                  panel.background = element_rect(fill = "white", colour = "grey90"),
                  strip.text.x = element_text(size = 12, colour = 'black'),
                  legend.key.width = unit(0.5, "cm"),
                  legend.key.height = unit(0.5, "cm"),
                  axis.text = element_text(color = "grey40"),
                  axis.ticks = element_blank(),
                  strip.text = element_text(colour = "black", size = 12),
                  panel.grid = element_line(color = "grey90"),
                  panel.grid.major = element_blank(),
                  plot.margin = margin(l = 5, b = 5, r = 5, t = 5),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold"),
                  strip.background = element_rect(colour = "grey90", fill ="white"),
                  legend.position = "right", legend.title = element_text(size = 12),
                  legend.text = element_text(size = 12),
                  legend.key = element_rect(fill = NA))

# colours
col1 <- "#AE445A"
col2 <- "#4B4376"
col3 <- "#B1C29E"
col4 <- "#FFE6A9"
col5 <- "#DEAA79"
col6 <- "#447D9B"
col7 <- "#9FC87E"
col8 <- "#FF9898"
col9 <- "#9B7EBD"
col10 <- "#309898"
