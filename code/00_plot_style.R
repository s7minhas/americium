# style publication figures

paper_plot_theme = function(base_size = 11) {
	ggplot2::theme_bw(base_size = base_size) +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			axis.line = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			legend.position = "top",
			plot.title = ggplot2::element_text(hjust = 0),
			plot.subtitle = ggplot2::element_text(hjust = 0),
			plot.caption = ggplot2::element_text(hjust = 0),
			strip.background = ggplot2::element_rect(
				fill = "black",
				color = "black"
			),
			strip.text.x = ggplot2::element_text(color = "white", hjust = 0),
			strip.text.y = ggplot2::element_text(color = "white", hjust = 0)
		)
}
