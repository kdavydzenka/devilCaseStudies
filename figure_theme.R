
library(ggplot2)
MY_THEME = ggplot2::theme(
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "gray90"),
  legend.position = "bottom",
  legend.title = element_text(face = "bold"),
  panel.grid.minor = element_blank()
)
