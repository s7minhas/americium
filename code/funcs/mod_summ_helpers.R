# prepare coefficient intervals and plots

coefplot_colors = c(
	'Positive' = rgb(54, 144, 192, maxColorValue = 255),
	'Negative' = rgb(222, 45, 38, maxColorValue = 255),
	'Positive at 90' = rgb(158, 202, 225, maxColorValue = 255),
	'Negative at 90' = rgb(252, 146, 114, maxColorValue = 255),
	'Not Significant' = rgb(150, 150, 150, maxColorValue = 255)
)

get_conf_ints = function(est_data, beta = 'Estimate', serror = 'Std. Error'){
	est_data$estimate = est_data[, beta]
	est_data$qt95_hi = est_data[, beta] + qnorm(0.975) * est_data[, serror]
	est_data$qt95_lo = est_data[, beta] - qnorm(0.975) * est_data[, serror]
	est_data$qt90_hi = est_data[, beta] + qnorm(0.95) * est_data[, serror]
	est_data$qt90_lo = est_data[, beta] - qnorm(0.95) * est_data[, serror]
	est_data
}

get_sig_levels = function(est_data){
	est_data$sig = 'Not Significant'
	est_data$sig[est_data$qt90_lo > 0 & est_data$qt95_lo < 0] = 'Positive at 90'
	est_data$sig[est_data$qt95_lo > 0] = 'Positive'
	est_data$sig[est_data$qt90_hi < 0 & est_data$qt95_hi > 0] = 'Negative at 90'
	est_data$sig[est_data$qt95_hi < 0] = 'Negative'
	est_data
}

coef_process = function(est_data, beta = 'Estimate', serror = 'Std. Error'){
	get_sig_levels(
		get_conf_ints(
			est_data = est_data,
			beta = beta,
			serror = serror
		)
	)
}

coef_viz = function(gg_data){
	ggplot(
		data = gg_data,
		aes(x = iv_clean, y = estimate, color = sig)
	) +
		geom_point() +
		geom_linerange(
			aes(ymin = qt95_lo, ymax = qt95_hi),
			alpha = 1,
			size = 0.7
		) +
		geom_linerange(
			aes(ymin = qt90_lo, ymax = qt90_hi),
			alpha = 1,
			size = 1
		) +
		scale_color_manual(values = coefplot_colors) +
		geom_hline(
			aes(yintercept = 0),
			color = 'grey50',
			linetype = 'dashed'
		) +
		labs(
			x = '',
			y = '',
			color = ''
		) +
		coord_flip() +
		facet_wrap(~mod_lab, nrow = 1, scales = 'free_x') +
		theme_bw() +
		theme(
			legend.position = 'top',
			axis.ticks = element_blank(),
			panel.border = element_blank(),
			strip.background = element_rect(fill = 'black', color = 'black'),
			strip.text.x = element_text(color = 'white', hjust = 0),
			strip.text.y = element_text(color = 'white', hjust = 0)
		)
}
