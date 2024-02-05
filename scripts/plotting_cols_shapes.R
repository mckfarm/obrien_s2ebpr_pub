# plotting params

library(ggplot2)

theme_black_box <- list(
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = NA))
)

reactor_cols <- c("#06d6a0", "#118ab2", "#073b4c")
carbon_cols <- c("#3a86ff", "#ff006e", "#8338ec")
mode_cols <- c("#ff9770", "#ff70a6")

scale_color_reactor <- scale_color_manual(values = reactor_cols, name = "Reactor")
scale_shape_reactor <- scale_shape_manual(values = c(21, 22, 23), name = "Reactor")

scale_color_reactor <- scale_color_manual(values = reactor_cols, name = "Reactor")
scale_shape_op_mode <- scale_shape_manual(values = c(16, 17), name = "Operation mode")
scale_color_perf <- scale_color_manual(values = carbon_cols, name = "Carbon")


scale_shape_carb <- scale_shape_manual(values = c(21, 19), name = "Carbon dose")
