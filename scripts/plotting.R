# plotting params

library(ggplot2)

reactor_cols <- c("#c44900", "#183a37", "#693A52")

scale_color_reactor <- scale_color_manual(values = reactor_cols, name = "Reactor")
scale_shape_reactor <- scale_shape_manual(values = c(21, 22, 23), name = "Reactor")

phases <- data.frame(x0 = ymd("2022-10-26"),
                     x1 = ymd("2023-01-24"),
                     x2 = ymd("2023-03-06"),
                     x3 = ymd("2023-03-22"),
                     x4 = ymd("2023-05-12"))

sbr3_lines <- list(
  geom_vline(xintercept = phases$x2, color="navyblue", linetype="dashed"),
  geom_vline(xintercept = phases$x3, color="navyblue", linetype="dashed")
)

x_axis_date <- scale_x_date(
  breaks="3 months",
  date_labels="%b-%y"
)

x_axis_date_short <- scale_x_date(
  breaks="1 months",
  date_labels="%b"
)


default_timeseries <- list(theme_bw(), 
                           labs(y = "Relative abundance [%]", x = "Date"), 
                           x_axis_date
)

default_timeseries_short <- list(theme_bw(), 
                           labs(y = "Relative abundance [%]", x = "Date"), 
                           x_axis_date_short
)
