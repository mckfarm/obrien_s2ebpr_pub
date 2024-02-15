phases <- data.frame(x0 = ymd("2022-10-26"),
                     x1 = ymd("2023-01-24"),
                     x2 = ymd("2023-03-06"),
                     x3 = ymd("2023-03-22"),
                     x4 = ymd("2023-05-12"),
                     red_C = ymd("2022-12-2"),
                     no_C = ymd("2023-05-12"))

# sbr3_lines <- list(
#   geom_vline(xintercept = phases$x2, color="navyblue", linetype="dashed"),
#   geom_vline(xintercept = phases$x3, color="navyblue", linetype="dashed")
# )
# 
# carbon_lines <- list(
#   geom_vline(xintercept = phases$x0, color="navyblue", linetype="dashed"),
#   geom_vline(xintercept = phases$red_C, color="navyblue", linetype="dashed"),
#   geom_vline(xintercept = phases$no_C, color="navyblue", linetype="dashed")
# )

x_axis_date <- scale_x_date(
  breaks="3 months",
  date_labels="%m/%y"
)

x_axis_date_short <- scale_x_date(
  breaks="2 months",
  date_labels="%m/%y"
)
