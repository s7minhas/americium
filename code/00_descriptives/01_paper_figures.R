# generate paper figures

rm(list = ls())
setup_candidates = c(
	file.path(here::here(), 'setup.R'),
	file.path(here::here(), 'code', 'setup.R'),
	file.path(getwd(), 'setup.R'),
	file.path(getwd(), 'code', 'setup.R'),
	file.path(dirname(getwd()), 'code', 'setup.R')
)
setup_candidates = setup_candidates[file.exists(setup_candidates)]
if(length(setup_candidates) == 0){ stop('Could not locate setup.R.') }
source(setup_candidates[1])
source(paste0(path_funcs, 'project_metadata.R'))

load_pkg(c('ggplot2', 'reshape2', 'ggrepel'))

# shared theme
theme_paper = theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.ticks = element_blank(),
		legend.position = 'top',
		strip.background = element_rect(fill = 'black', color = 'black'),
		strip.text = element_text(color = 'white', hjust = 0)
	)

save_fig = function(plot, name, width, height){
	ggsave(paste0(path_graphics, name, '.pdf'), plot, width = width, height = height)
	ggsave(paste0(path_graphics, name, '.png'), plot, width = width, height = height, dpi = 600)
}

type_levels = amer_category_labels
type_colors = amer_category_colors
history_type_levels = c(type_levels, 'Other States')
history_type_colors = c(type_colors, 'Other States' = 'grey78')

# cooperation space
load(paste0(path_out, 'snap_fit_un_voting.rda'))
u_arr = fit$U
u_sd = fit$U_sd
label_focus = c('USA','DEU','JPN','RUS','CHN','IND','TUR','BRA')

observed_actors = function(yr){
	y_yr = fit$Y[, , yr]
	rownames(y_yr)[rowSums(!is.na(y_yr)) + colSums(!is.na(y_yr)) > 0]
}

historical_actor_type = function(codes){
	out = amer_actor_type(codes)
	out[codes == 'SUN'] = 'Strategic Rivals'
	out[codes == 'GER'] = 'Core Allies'
	out
}

fit_alignment = function(target, reference, common_actors){
	ref = reference[common_actors, , drop = FALSE]
	tgt = target[common_actors, , drop = FALSE]
	ref_center = colMeans(ref)
	tgt_center = colMeans(tgt)
	ref_centered = sweep(ref, 2, ref_center)
	tgt_centered = sweep(tgt, 2, tgt_center)
	rot = svd(t(tgt_centered) %*% ref_centered)
	q = rot$u %*% t(rot$v)
	aligned = sweep(target, 2, tgt_center)
	aligned = aligned %*% q
	list(
		aligned = sweep(aligned, 2, ref_center, '+'),
		rotation = q,
		target_center = tgt_center,
		reference_center = ref_center
	)
}

align_to_reference = function(target, reference, common_actors){
	fit_alignment(target, reference, common_actors)$aligned
}

latent_distance = function(a, b){
	sqrt(sum((a - b)^2))
}

mean_pair_distance = function(pos, actors_a, actors_b = NULL){
	actors_a = intersect(actors_a, rownames(pos))
	if(is.null(actors_b)){
		if(length(actors_a) < 2){ return(NA_real_) }
		pairs = combn(actors_a, 2)
		return(mean(apply(pairs, 2, function(z){
			latent_distance(pos[z[1], ], pos[z[2], ])
		}), na.rm = TRUE))
	}
	actors_b = intersect(actors_b, rownames(pos))
	if(length(actors_a) == 0 || length(actors_b) == 0){ return(NA_real_) }
	mean(as.vector(sapply(actors_a, function(i){
		sapply(actors_b, function(j){
			latent_distance(pos[i, ], pos[j, ])
		})
	})), na.rm = TRUE)
}

# align 2025 to 2024
pos_24 = u_arr[, , '2024']
pos_25_raw = u_arr[, , '2025']
common_24_25 = intersect(observed_actors('2024'), observed_actors('2025'))
space_alignment = fit_alignment(pos_25_raw, pos_24, common_24_25)
pos_25 = space_alignment$aligned
q_rot = space_alignment$rotation

build_year = function(pos, yr){
	df = data.frame(actor = rownames(pos), x = pos[, 1], y = pos[, 2], year = yr)
	if(!is.null(u_sd)){
		sd_yr = u_sd[, , yr]
		# rotate 2025 sds
		if(yr == '2025'){ sd_yr = sqrt((sd_yr^2) %*% (q_rot^2)) }
		df$sd_x = sd_yr[, 1]
		df$sd_y = sd_yr[, 2]
	} else { df$sd_x = NA; df$sd_y = NA }
	df
}
space_df = rbind(build_year(pos_24, '2024'), build_year(pos_25, '2025'))
space_df$type = amer_actor_type(space_df$actor)
space_df = space_df[!is.na(space_df$type), ]
space_df$type = factor(space_df$type, levels = type_levels)
space_df$label = ifelse(space_df$actor %in% label_focus, space_df$actor, NA)

# build uncertainty ellipses
fig_space = ggplot(space_df, aes(x = x, y = y, color = type))
if(!is.null(u_sd)){
	ang = seq(0, 2*pi, length.out = 40)
	ell = do.call(rbind, lapply(seq_len(nrow(space_df)), function(k){
		r = space_df[k, ]
		data.frame(
			x = r$x + 2 * r$sd_x * cos(ang),
			y = r$y + 2 * r$sd_y * sin(ang),
			grp = paste(r$actor, r$year),
			type = r$type, year = r$year
		)
	}))
	fig_space = fig_space +
		geom_polygon(data = ell, aes(x = x, y = y, group = grp, color = type),
			fill = NA, alpha = 0.4, linewidth = 0.2, show.legend = FALSE)
}
fig_space = fig_space +
	geom_point(alpha = 0.7, size = 1.2) +
	geom_text_repel(
		data = space_df[space_df$actor %in% label_focus, ],
		aes(label = label),
		size = 2.8, max.overlaps = Inf, box.padding = 0.25,
		point.padding = 0.2, min.segment.length = 0, segment.alpha = 0.35,
		seed = seed, show.legend = FALSE, na.rm = TRUE
	) +
	facet_wrap(~ year, ncol = 2) +
	scale_color_manual(values = type_colors, name = NULL) +
	labs(x = 'Latent Dimension 1', y = 'Latent Dimension 2') +
	theme_paper
save_fig(fig_space, 'cooperation_space', 9, 5)

# historical space shifts
history_events = data.frame(
	episode = c(
		'Post-Cold War\n1989 to 1992',
		'Post-9/11 and Bush\n2000 to 2002',
		'Trump I\n2016 to 2020',
		'2025 rupture\n2024 to 2025'
	),
	before = c('1989','2000','2016','2024'),
	after = c('1992','2002','2020','2025'),
	stringsAsFactors = FALSE
)
history_events$episode = factor(history_events$episode, levels = history_events$episode)
history_focus = unique(c(amer_focal_actors, 'SUN', 'GER'))
history_labels = c('USA','DEU','GER','SUN','RUS','CHN','JPN','IND','TUR','BRA')

history_parts = lapply(seq_len(nrow(history_events)), function(k){
	ev = history_events[k, ]
	pos_before = u_arr[, , ev$before]
	pos_after = u_arr[, , ev$after]
	obs_before = observed_actors(ev$before)
	obs_after = observed_actors(ev$after)
	common_obs = intersect(obs_before, obs_after)
	pos_after_aligned = align_to_reference(pos_after, pos_before, common_obs)

	point_before = data.frame(
		episode = ev$episode,
		actor = obs_before,
		stage = 'Before',
		year = ev$before,
		x = pos_before[obs_before, 1],
		y = pos_before[obs_before, 2],
		stringsAsFactors = FALSE
	)
	point_after = data.frame(
		episode = ev$episode,
		actor = obs_after,
		stage = 'After',
		year = ev$after,
		x = pos_after_aligned[obs_after, 1],
		y = pos_after_aligned[obs_after, 2],
		stringsAsFactors = FALSE
	)
	points = rbind(point_before, point_after)
	points$type = historical_actor_type(points$actor)
	points$type_plot = ifelse(is.na(points$type), 'Other States', points$type)
	points$type_plot = factor(points$type_plot, levels = history_type_levels)
	points$stage = factor(points$stage, levels = c('Before','After'))
	points$label = ifelse(points$actor %in% history_labels, points$actor, NA)

	seg_actors = intersect(history_focus, common_obs)
	segments = data.frame(
		episode = ev$episode,
		actor = seg_actors,
		x = pos_before[seg_actors, 1],
		y = pos_before[seg_actors, 2],
		xend = pos_after_aligned[seg_actors, 1],
		yend = pos_after_aligned[seg_actors, 2],
		stringsAsFactors = FALSE
	)
	segments$type = historical_actor_type(segments$actor)
	segments = segments[!is.na(segments$type), ]
	segments$type = factor(segments$type, levels = type_levels)

	list(points = points, segments = segments)
})

history_points = do.call(rbind, lapply(history_parts, `[[`, 'points'))
history_segments = do.call(rbind, lapply(history_parts, `[[`, 'segments'))
history_points_plot = do.call(rbind, lapply(split(history_points, history_points$episode), function(df){
	x_keep = quantile(df$x, c(0.02, 0.98), na.rm = TRUE)
	y_keep = quantile(df$y, c(0.02, 0.98), na.rm = TRUE)
	keep = df$type_plot != 'Other States' |
		(df$x >= x_keep[1] & df$x <= x_keep[2] & df$y >= y_keep[1] & df$y <= y_keep[2])
	df[keep, ]
}))
history_labels_df = history_points[
	history_points$actor %in% history_labels &
		(history_points$stage == 'After' | history_points$actor %in% c('SUN','GER')),
]
history_labels_df = history_labels_df[history_labels_df$type_plot != 'Other States', ]

fig_history = ggplot() +
	geom_point(
		data = history_points_plot[history_points_plot$type_plot == 'Other States', ],
		aes(x = x, y = y, shape = stage),
		color = 'grey75', size = 0.65, alpha = 0.35
	) +
	geom_segment(
		data = history_segments,
		aes(x = x, y = y, xend = xend, yend = yend, color = type),
		arrow = grid::arrow(length = grid::unit(0.055, 'inches')),
		linewidth = 0.35, alpha = 0.65
	) +
	geom_point(
		data = history_points_plot[history_points_plot$type_plot != 'Other States', ],
		aes(x = x, y = y, color = type_plot, shape = stage),
		size = 1.15, alpha = 0.9
	) +
	geom_text_repel(
		data = history_labels_df,
		aes(x = x, y = y, label = label, color = type_plot),
		size = 2.3, max.overlaps = Inf, box.padding = 0.2,
		point.padding = 0.15, min.segment.length = 0, segment.alpha = 0.3,
		seed = seed, show.legend = FALSE, na.rm = TRUE
	) +
	facet_wrap(~ episode, ncol = 2, scales = 'free') +
	scale_color_manual(values = history_type_colors, name = 'Actor category') +
	scale_shape_manual(values = c('Before' = 1, 'After' = 16), name = NULL) +
	labs(x = 'Latent Dimension 1', y = 'Latent Dimension 2') +
	theme_paper +
	theme(
		legend.position = 'top',
		strip.text = element_text(color = 'white', hjust = 0.5),
		axis.text = element_text(size = 7)
	)
save_fig(fig_history, 'historical_space_shifts', 8.5, 6.2)

# compact historical shifts
centroid_groups = list(
	'Hegemon' = 'USA',
	'Core Allies' = unique(c(amer_actor_cats$allies, 'GER')),
	'Strategic Rivals' = unique(c(amer_actor_cats$rivals, 'SUN')),
	'Pivotal Hedgers' = amer_actor_cats$hedgers
)
centroid_events = data.frame(
	episode = c(
		'Soviet collapse\n1989 to 1992',
		'9/11 and Bush\n2000 to 2002',
		'Trump I\n2016 to 2020',
		'2025 rupture\n2024 to 2025'
	),
	before = c('1989','2000','2016','2024'),
	after = c('1992','2002','2020','2025'),
	stringsAsFactors = FALSE
)
centroid_events$episode = factor(centroid_events$episode, levels = centroid_events$episode)

centroid_parts = lapply(seq_len(nrow(centroid_events)), function(k){
	ev = centroid_events[k, ]
	pos_before = u_arr[, , ev$before]
	pos_after = u_arr[, , ev$after]
	obs_before = observed_actors(ev$before)
	obs_after = observed_actors(ev$after)
	common_obs = intersect(obs_before, obs_after)
	pos_after_aligned = align_to_reference(pos_after, pos_before, common_obs)

	points = do.call(rbind, lapply(names(centroid_groups), function(grp){
		actors_before = intersect(centroid_groups[[grp]], obs_before)
		actors_after = intersect(centroid_groups[[grp]], obs_after)
		if(length(actors_before) == 0 || length(actors_after) == 0){ return(NULL) }
		data.frame(
			episode = ev$episode,
			group = grp,
			stage = c('Before', 'After'),
			year = c(ev$before, ev$after),
			x = c(mean(pos_before[actors_before, 1]), mean(pos_after_aligned[actors_after, 1])),
			y = c(mean(pos_before[actors_before, 2]), mean(pos_after_aligned[actors_after, 2])),
			n = c(length(actors_before), length(actors_after)),
			stringsAsFactors = FALSE
		)
	}))
	segments = do.call(rbind, lapply(split(points, points$group), function(df){
		if(nrow(df) != 2){ return(NULL) }
		before = df[df$stage == 'Before', ]
		after = df[df$stage == 'After', ]
		data.frame(
			episode = ev$episode,
			group = before$group,
			x = before$x, y = before$y,
			xend = after$x, yend = after$y,
			stringsAsFactors = FALSE
		)
	}))
	background = data.frame(
		episode = ev$episode,
		x = pos_after_aligned[obs_after, 1],
		y = pos_after_aligned[obs_after, 2],
		stringsAsFactors = FALSE
	)
	list(points = points, segments = segments, background = background)
})

centroid_points = do.call(rbind, lapply(centroid_parts, `[[`, 'points'))
centroid_segments = do.call(rbind, lapply(centroid_parts, `[[`, 'segments'))
centroid_background = do.call(rbind, lapply(centroid_parts, `[[`, 'background'))
centroid_points$episode = factor(centroid_points$episode, levels = levels(centroid_events$episode))
centroid_segments$episode = factor(centroid_segments$episode, levels = levels(centroid_events$episode))
centroid_background$episode = factor(centroid_background$episode, levels = levels(centroid_events$episode))
centroid_points$group = factor(centroid_points$group, levels = type_levels)
centroid_segments$group = factor(centroid_segments$group, levels = type_levels)
centroid_points$stage = factor(centroid_points$stage, levels = c('Before', 'After'))
centroid_labels = centroid_points[centroid_points$stage == 'After', ]
centroid_labels$label = c(
	'Hegemon' = 'U.S.',
	'Core Allies' = 'Core allies',
	'Strategic Rivals' = 'Rivals',
	'Pivotal Hedgers' = 'Hedgers'
)[as.character(centroid_labels$group)]

fig_centroids = ggplot() +
	geom_point(
		data = centroid_background,
		aes(x = x, y = y),
		color = 'grey78', size = 0.55, alpha = 0.3
	) +
	geom_segment(
		data = centroid_segments,
		aes(x = x, y = y, xend = xend, yend = yend, color = group),
		arrow = grid::arrow(length = grid::unit(0.07, 'inches')),
		linewidth = 0.6, alpha = 0.75
	) +
	geom_point(
		data = centroid_points,
		aes(x = x, y = y, color = group, shape = stage),
		size = 2.7, stroke = 0.5
	) +
	geom_text_repel(
		data = centroid_labels,
		aes(x = x, y = y, label = label, color = group),
		size = 2.55, max.overlaps = Inf, box.padding = 0.25,
		point.padding = 0.18, min.segment.length = 0, segment.alpha = 0.35,
		seed = seed, show.legend = FALSE
	) +
	facet_wrap(~ episode, ncol = 2, scales = 'free') +
	scale_color_manual(values = type_colors, name = 'Group') +
	scale_shape_manual(values = c('Before' = 1, 'After' = 16), name = NULL) +
	labs(x = 'Latent Dimension 1', y = 'Latent Dimension 2') +
	theme_paper +
	theme(
		legend.position = 'top',
		strip.text = element_text(color = 'white', hjust = 0.5),
		axis.text = element_text(size = 7)
	)
save_fig(fig_centroids, 'historical_centroid_shifts', 8.5, 5.6)

# allied bridge distances
bridge_europe = c('FRA','DEU','GER','GBR','ITA','NLD','ESP','BEL','NOR','DNK','POL')
bridge_indo = c('JPN','KOR','AUS')

bridge_metrics = do.call(rbind, lapply(dimnames(u_arr)[[3]], function(yr){
	pos = u_arr[, , yr]
	obs = observed_actors(yr)
	europe = intersect(bridge_europe, obs)
	indo = intersect(bridge_indo, obs)
	data.frame(
		year = as.integer(yr),
		us_to_europe = mean_pair_distance(pos, 'USA', europe),
		us_to_indo = mean_pair_distance(pos, 'USA', indo),
		europe_to_indo = mean_pair_distance(pos, europe, indo),
		japan_to_europe = mean_pair_distance(pos, 'JPN', europe),
		korea_to_europe = mean_pair_distance(pos, 'KOR', europe),
		australia_to_europe = mean_pair_distance(pos, 'AUS', europe),
		stringsAsFactors = FALSE
	)
}))

bridge_long = reshape2::melt(bridge_metrics, id.vars = 'year',
	variable.name = 'metric', value.name = 'distance')
bridge_labels = c(
	us_to_europe = 'U.S. to Europe',
	us_to_indo = 'U.S. to Japan/Korea/Australia',
	europe_to_indo = 'Europe to Japan/Korea/Australia',
	japan_to_europe = 'Japan to Europe',
	korea_to_europe = 'South Korea to Europe',
	australia_to_europe = 'Australia to Europe'
)
bridge_long$series = bridge_labels[as.character(bridge_long$metric)]
bridge_long$series = factor(bridge_long$series, levels = unname(bridge_labels))
bridge_long$panel = ifelse(
	bridge_long$metric %in% c('us_to_europe','us_to_indo','europe_to_indo'),
	'U.S. distance to allied wings',
	'Indo-Pacific allies to Europe'
)
bridge_long$panel = factor(bridge_long$panel,
	levels = c('U.S. distance to allied wings', 'Indo-Pacific allies to Europe'))
bridge_colors = c(
	'U.S. to Europe' = '#d62728',
	'U.S. to Japan/Korea/Australia' = '#b44e19',
	'Europe to Japan/Korea/Australia' = '#1f77b4',
	'Japan to Europe' = '#1f77b4',
	'South Korea to Europe' = '#17becf',
	'Australia to Europe' = 'grey40'
)
bridge_end = bridge_long[bridge_long$year == 2025, ]

fig_bridge = ggplot(bridge_long, aes(x = year, y = distance, color = series)) +
	geom_vline(xintercept = 2025, linetype = 'dashed', color = 'grey65', linewidth = 0.3) +
	geom_line(linewidth = 0.65) +
	geom_point(data = bridge_end, size = 1.8, show.legend = FALSE) +
	geom_text_repel(
		data = bridge_end, aes(label = series),
		nudge_x = 1.3, direction = 'y', hjust = 0, size = 2.55,
		min.segment.length = 0, segment.alpha = 0.35, box.padding = 0.2,
		seed = seed, show.legend = FALSE
	) +
	facet_wrap(~ panel, ncol = 1, scales = 'free_y') +
	scale_color_manual(values = bridge_colors, guide = 'none') +
	scale_x_continuous(breaks = seq(1990, 2025, 5)) +
	labs(x = 'Year', y = 'Mean Latent Distance') +
	coord_cartesian(xlim = c(1989, 2028), clip = 'off') +
	theme_paper +
	theme(plot.margin = margin(5.5, 70, 5.5, 5.5))
save_fig(fig_bridge, 'allied_bridge_distances', 7.5, 5.4)

# allied bridge schematic
bridge_nodes = data.frame(
	node = c('U.S.', 'Europe', 'JPN/KOR/AUS'),
	x = c(0.5, 0.18, 0.82),
	y = c(0.86, 0.18, 0.18),
	color = c('#d62728', '#1f77b4', '#2ca02c'),
	stringsAsFactors = FALSE
)
bridge_edge_lookup = data.frame(
	metric = c('us_to_europe', 'us_to_indo', 'europe_to_indo'),
	relation = c('U.S.-Europe', 'U.S.-JPN/KOR/AUS', 'Europe-JPN/KOR/AUS'),
	from = c('U.S.', 'U.S.', 'Europe'),
	to = c('Europe', 'JPN/KOR/AUS', 'JPN/KOR/AUS'),
	color = c('#d62728', '#b44e19', '#1f77b4'),
	stringsAsFactors = FALSE
)
bridge_schematic_edges = merge(
	bridge_long[bridge_long$year %in% c(2024, 2025) &
		bridge_long$metric %in% bridge_edge_lookup$metric, c('year','metric','distance')],
	bridge_edge_lookup,
	by = 'metric'
)
bridge_schematic_edges = merge(
	bridge_schematic_edges, bridge_nodes[, c('node','x','y')],
	by.x = 'from', by.y = 'node'
)
bridge_schematic_edges = merge(
	bridge_schematic_edges, bridge_nodes[, c('node','x','y')],
	by.x = 'to', by.y = 'node', suffixes = c('', 'end')
)
bridge_schematic_edges$year = factor(bridge_schematic_edges$year, levels = c(2024, 2025))
bridge_schematic_edges$label = sprintf('%.3f', bridge_schematic_edges$distance)
bridge_schematic_edges$label_x = (bridge_schematic_edges$x + bridge_schematic_edges$xend) / 2
bridge_schematic_edges$label_y = (bridge_schematic_edges$y + bridge_schematic_edges$yend) / 2
bridge_schematic_edges$label_y = bridge_schematic_edges$label_y +
	ifelse(bridge_schematic_edges$relation == 'Europe-JPN/KOR/AUS', -0.07, 0.04)
bridge_schematic_edges$relation = factor(bridge_schematic_edges$relation,
	levels = bridge_edge_lookup$relation)
bridge_schematic_nodes = do.call(rbind, lapply(levels(bridge_schematic_edges$year), function(yr){
	df = bridge_nodes
	df$year = yr
	df
}))
bridge_schematic_nodes$year = factor(bridge_schematic_nodes$year,
	levels = levels(bridge_schematic_edges$year))

fig_bridge_space = ggplot() +
	geom_segment(
		data = bridge_schematic_edges,
		aes(x = x, y = y, xend = xend, yend = yend,
			linewidth = distance, color = relation),
		alpha = 0.75, lineend = 'round'
	) +
	geom_text(
		data = bridge_schematic_edges,
		aes(x = label_x, y = label_y, label = label, color = relation),
		size = 3.1, show.legend = FALSE
	) +
	geom_point(
		data = bridge_schematic_nodes,
		aes(x = x, y = y, fill = node),
		shape = 21, color = 'white', stroke = 0.45, size = 5
	) +
	geom_text_repel(
		data = bridge_schematic_nodes,
		aes(x = x, y = y, label = node),
		size = 3.2, fontface = 'bold', color = 'black',
		box.padding = 0.3, point.padding = 0.22, min.segment.length = 0,
		segment.alpha = 0.25, seed = seed + 1, show.legend = FALSE
	) +
	facet_wrap(~ year, ncol = 2) +
	scale_color_manual(values = setNames(bridge_edge_lookup$color, bridge_edge_lookup$relation),
		guide = 'none') +
	scale_fill_manual(values = setNames(bridge_nodes$color, bridge_nodes$node),
		guide = 'none') +
	scale_linewidth(range = c(0.6, 3.3), guide = 'none') +
	labs(x = NULL, y = NULL) +
	coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = 'off') +
	theme_void() +
	theme(
		plot.background = element_rect(fill = 'white', color = NA),
		panel.background = element_rect(fill = 'white', color = NA),
		strip.background = element_rect(fill = 'black', color = 'black'),
		strip.text = element_text(color = 'white', size = 11, hjust = 0.5),
		panel.spacing.x = grid::unit(0.45, 'lines'),
		plot.margin = margin(5.5, 16, 5.5, 16)
	)
save_fig(fig_bridge_space, 'allied_bridge_space', 7.2, 3.7)

# 2025 snaps by category
sp_2025 = fit$snap_prob[, '2025']
snap_df = data.frame(actor = names(sp_2025), snap = as.numeric(sp_2025))
snap_df$type = amer_actor_type(snap_df$actor)
snap_df = snap_df[!is.na(snap_df$type), ]
snap_df$type = factor(snap_df$type, levels = type_levels)
snap_df = snap_df[order(snap_df$type, snap_df$snap, snap_df$actor), ]
snap_df$observed_2025 = snap_df$actor %in% observed_actors('2025')
snap_band_levels = c('No 2025 UN ties', '< 0.50', '0.50-0.89', '0.90-0.98', '0.99-1.00')
snap_df$band = cut(
	snap_df$snap,
	breaks = c(-Inf, 0.5, 0.9, 0.99, Inf),
	labels = snap_band_levels[-1],
	right = FALSE
)
snap_df$band = as.character(snap_df$band)
snap_df$band[!snap_df$observed_2025] = 'No 2025 UN ties'
snap_df$band = factor(snap_df$band, levels = snap_band_levels)

snap_counts = as.data.frame(table(snap_df$type, snap_df$band), stringsAsFactors = FALSE)
names(snap_counts) = c('type', 'band', 'n')
snap_counts$type = factor(snap_counts$type, levels = type_levels)
snap_counts$band = factor(snap_counts$band, levels = snap_band_levels)
snap_counts$n = as.numeric(snap_counts$n)
snap_counts$label_color = ifelse(snap_counts$band == '0.99-1.00', 'white', 'grey20')

snap_summary = aggregate(
	cbind(observed_2025 = observed_2025, above = observed_2025 & snap >= 0.5) ~ type,
	snap_df,
	sum
)
snap_summary$total = as.numeric(table(snap_df$type)[as.character(snap_summary$type)])
snap_summary$mean_observed = vapply(as.character(snap_summary$type), function(tt){
	mean(snap_df$snap[snap_df$type == tt & snap_df$observed_2025], na.rm = TRUE)
}, numeric(1))
snap_summary$label = paste0(
	snap_summary$above, '/', snap_summary$observed_2025,
	' observed > 0.5\nmean = ', sprintf('%.2f', snap_summary$mean_observed)
)
snap_summary$y = snap_summary$total + 0.65

snap_band_colors = c(
	'No 2025 UN ties' = 'grey82',
	'< 0.50' = 'grey95',
	'0.50-0.89' = '#fdae61',
	'0.90-0.98' = '#74a9cf',
	'0.99-1.00' = '#2b8cbe'
)

fig_snap = ggplot(snap_counts, aes(x = type, y = n, fill = band)) +
	geom_col(width = 0.68, color = 'white', linewidth = 0.35) +
	geom_text(
		data = snap_counts[snap_counts$n > 0, ],
		aes(label = n, color = label_color),
		position = position_stack(vjust = 0.5),
		size = 3.1, fontface = 'bold', show.legend = FALSE
	) +
	geom_text(
		data = snap_summary,
		aes(x = type, y = y, label = label),
		inherit.aes = FALSE,
		size = 2.7, color = 'grey20', lineheight = 0.9
	) +
	scale_fill_manual(values = snap_band_colors, drop = FALSE, name = 'Posterior snap band') +
	scale_color_identity() +
	scale_y_continuous(
		limits = c(0, max(snap_summary$y) + 0.9),
		breaks = seq(0, 14, 2),
		expand = c(0, 0)
	) +
	labs(x = NULL, y = 'Number of Focal Actors') +
	theme_paper +
	theme(
		legend.position = 'top',
		axis.text.x = element_text(size = 8.5),
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		panel.grid.minor.y = element_blank(),
		plot.margin = margin(5.5, 14, 5.5, 14)
	)
save_fig(fig_snap, 'snap_by_state', 7.4, 4.2)

# focal snap heatmap
snap_years = setdiff(colnames(fit$snap_prob), colnames(fit$snap_prob)[1])
snap_actor_order = intersect(unique(unlist(amer_actor_cats, use.names = FALSE)),
	rownames(fit$snap_prob))
snap_hist = as.data.frame(as.table(fit$snap_prob[snap_actor_order, snap_years, drop = FALSE]),
	stringsAsFactors = FALSE)
names(snap_hist) = c('actor', 'year', 'snap')
snap_hist$year = as.integer(as.character(snap_hist$year))
snap_hist$type = factor(amer_actor_type(snap_hist$actor), levels = type_levels)
snap_hist$actor_plot = factor(snap_hist$actor, levels = rev(snap_actor_order))

snap_heat_events = data.frame(year = c(1991, 2001, 2008, 2014, 2016, 2022, 2025))
fig_snap_hist = ggplot(snap_hist, aes(x = year, y = actor_plot, fill = snap)) +
	geom_tile(width = 0.92, height = 0.88) +
	geom_vline(data = snap_heat_events, aes(xintercept = year),
		linetype = 'dashed', color = 'grey25', linewidth = 0.45, alpha = 0.75,
		inherit.aes = FALSE) +
	facet_grid(type ~ ., scales = 'free_y', space = 'free_y') +
	scale_fill_gradientn(
		colors = c('white', '#cfe8f3', '#74a9cf', '#fdae61', '#b2182b'),
		limits = c(0, 1), breaks = seq(0, 1, 0.25),
		name = 'Posterior\nSnap Probability'
	) +
	guides(fill = guide_colorbar(
		title.position = 'top', title.hjust = 0.5,
		barwidth = grid::unit(2.8, 'inches'),
		barheight = grid::unit(0.18, 'inches')
	)) +
	scale_x_continuous(breaks = c(seq(1990, 2020, 5), 2025), expand = c(0.01, 0.01)) +
	labs(x = 'Year', y = NULL) +
	theme_paper +
	theme(
		legend.position = 'top',
		legend.title = element_text(size = 8),
		legend.text = element_text(size = 7),
		axis.text.y = element_text(size = 7),
		strip.background = element_blank(),
		strip.text.y = element_blank(),
		panel.spacing.y = grid::unit(0.35, 'lines')
	)
save_fig(fig_snap_hist, 'snap_history_heatmap', 7.5, 7.2)

# period snap summary
period_windows = data.frame(
	period = c(
		'Post-Cold War\n1990-1992',
		'Late 1990s\n1996-1998',
		'9/11 and Bush\n2001-2003',
		'Financial crisis\n2008-2009',
		'Crimea aftermath\n2014-2015',
		'Trump I\n2017-2020',
		'Ukraine invasion\n2022-2023',
		'2025 rupture\n2025'
	),
	years = I(list(
		1990:1992,
		1996:1998,
		2001:2003,
		2008:2009,
		2014:2015,
		2017:2020,
		2022:2023,
		2025
	)),
	stringsAsFactors = FALSE
)
period_windows$period = factor(period_windows$period, levels = period_windows$period)

period_cat = do.call(rbind, lapply(seq_len(nrow(period_windows)), function(k){
	period = as.character(period_windows$period[k])
	yrs = intersect(as.character(period_windows$years[[k]]), colnames(fit$snap_prob))
	do.call(rbind, lapply(names(amer_actor_cats), function(category_key){
		actors = intersect(amer_actor_cats[[category_key]], rownames(fit$snap_prob))
		vals = vapply(yrs, function(yr){
			mean(fit$snap_prob[actors, yr], na.rm = TRUE)
		}, numeric(1))
		data.frame(
			period = period,
			category = amer_category_labels[[category_key]],
			max_snap = max(vals, na.rm = TRUE),
			max_year = names(vals)[which.max(vals)],
			stringsAsFactors = FALSE
		)
	}))
}))
period_cat$period = factor(period_cat$period, levels = levels(period_windows$period))
period_cat$category = factor(period_cat$category, levels = rev(type_levels))
period_cat$label = paste0(sprintf('%.2f', period_cat$max_snap), '\n', period_cat$max_year)
period_cat$label_color = ifelse(period_cat$max_snap > 0.62, 'white', 'black')

fig_periods = ggplot(period_cat, aes(x = period, y = category, fill = max_snap)) +
	geom_tile(color = 'white', linewidth = 0.7) +
	geom_text(aes(label = label, color = label_color),
		size = 2.6, lineheight = 0.9) +
	scale_fill_gradientn(
		colors = c('white', '#d7eef3', '#8bbbd7', '#fdae61', '#b2182b'),
		limits = c(0, 1), breaks = seq(0, 1, 0.25),
		name = 'Maximum Mean\nSnap Probability'
	) +
	guides(fill = guide_colorbar(
		title.position = 'top', title.hjust = 0.5,
		barwidth = grid::unit(2.6, 'inches'),
		barheight = grid::unit(0.18, 'inches')
	)) +
	scale_color_identity() +
	labs(x = NULL, y = NULL) +
	theme_paper +
	theme(
		legend.position = 'top',
		legend.title = element_text(size = 8),
		legend.text = element_text(size = 7),
		axis.text.x = element_text(size = 7.5, hjust = 0.5),
		axis.text.y = element_text(size = 8),
		panel.grid = element_blank()
	)
save_fig(fig_periods, 'snap_periods_by_category', 8.5, 4.2)

# transition diagnostic
rho_hat = mean(fit$rho_uv_aligned, na.rm = TRUE)
sigma_hat = mean(fit$sigma_uv, na.rm = TRUE)
transition_df = do.call(rbind, lapply(seq(2, length(dimnames(u_arr)[[3]])), function(k){
	yr = dimnames(u_arr)[[3]][k]
	yr_prev = dimnames(u_arr)[[3]][k - 1]
	pos = u_arr[, , yr]
	pos_prev = u_arr[, , yr_prev]
	actors = Reduce(intersect, list(rownames(pos), observed_actors(yr), observed_actors(yr_prev)))
	ar_resid = sqrt(rowSums((pos[actors, , drop = FALSE] -
		rho_hat * pos_prev[actors, , drop = FALSE])^2))
	movement = sqrt(rowSums((pos[actors, , drop = FALSE] -
		pos_prev[actors, , drop = FALSE])^2))
	data.frame(
		actor = actors,
		year = as.integer(yr),
		movement = movement,
		drift_units = ar_resid / sigma_hat,
		snap = fit$snap_prob[actors, yr],
		type = amer_actor_type(actors),
		stringsAsFactors = FALSE
	)
}))
transition_df = transition_df[is.finite(transition_df$drift_units) &
	is.finite(transition_df$snap), ]
transition_df$drift_units_plot = pmax(transition_df$drift_units, 0.05)
transition_df$type = factor(transition_df$type, levels = type_levels)
transition_focus = transition_df[
	transition_df$year == 2025 & transition_df$actor %in% amer_focal_actors,
]
transition_label_actors = c('USA','FRA','DEU','GBR','JPN','AUS','KOR','RUS','CHN','IND','TUR','BRA')
transition_labels = transition_focus[transition_focus$actor %in% transition_label_actors, ]

fig_transition = ggplot(transition_df, aes(x = drift_units_plot, y = snap)) +
	geom_point(color = 'grey68', alpha = 0.22, size = 0.55) +
	geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'grey45', linewidth = 0.3) +
	geom_vline(xintercept = 1, linetype = 'dotted', color = 'grey55', linewidth = 0.3) +
	geom_point(
		data = transition_focus,
		aes(color = type),
		size = 1.8, alpha = 0.9
	) +
	geom_text_repel(
		data = transition_labels,
		aes(label = actor, color = type),
		size = 2.55, max.overlaps = Inf, box.padding = 0.25,
		point.padding = 0.16, min.segment.length = 0, segment.alpha = 0.35,
		seed = seed, show.legend = FALSE
	) +
	scale_x_log10(
		breaks = c(0.1, 0.3, 1, 3, 10, 30, 100),
		labels = c('0.1', '0.3', '1', '3', '10', '30', '100')
	) +
	scale_y_continuous(breaks = seq(0, 1, 0.25)) +
	scale_color_manual(values = type_colors, name = '2025 focal actors') +
	labs(
		x = 'Distance from AR-Drift Prediction (Drift SDs)',
		y = 'Posterior Snap Probability'
	) +
	coord_cartesian(ylim = c(0, 1.03), clip = 'off') +
	theme_paper +
	theme(plot.margin = margin(5.5, 18, 5.5, 5.5))
save_fig(fig_transition, 'snap_transition_diagnostic', 7.2, 4.8)

# rupture index
load(paste0(path_out, 'rupture_index.rda'))
calib = data.frame(
	year = as.integer(names(rupture_index$un_voting$index_by_year)),
	index = as.numeric(rupture_index$un_voting$index_by_year)
)
calib = calib[is.finite(calib$index), ]
calib$highlight = ifelse(calib$year == 2025, '2025', 'Other Years')
event_df = data.frame(
	year = c(1991, 2001, 2008, 2014, 2016, 2022, 2025),
	label = c(
		'Soviet\ncollapse',
		'9/11',
		'Financial\ncrisis',
		'Crimea\nannexed',
		'Trump\nelected',
		'Russia invades\nUkraine',
		'2025 U.S.\nbreak'
	),
	y = c(0.96, 0.86, 0.76, 0.66, 0.56, 0.78, 0.91)
)
fig_calib = ggplot(calib, aes(x = year, y = index)) +
	geom_vline(data = event_df, aes(xintercept = year), inherit.aes = FALSE,
		linetype = 'dashed', color = 'grey65', linewidth = 0.3) +
	geom_line() +
	geom_point(aes(color = highlight), size = 1.6) +
	geom_text_repel(
		data = event_df, aes(x = year, y = y, label = label),
		inherit.aes = FALSE, size = 2.6, lineheight = 0.9,
		min.segment.length = 0, box.padding = 0.2, point.padding = 0.15,
		segment.color = 'grey55', segment.alpha = 0.7, seed = seed,
		max.overlaps = Inf
	) +
	scale_color_manual(values = c('2025' = '#d62728', 'Other Years' = 'black'), guide = 'none') +
	labs(x = 'Year', y = 'Focal-Actor Rupture Index') +
	coord_cartesian(xlim = c(1989, 2026), ylim = c(0, 1.05), clip = 'off') +
	theme_paper +
	theme(plot.margin = margin(5.5, 16, 5.5, 5.5))
save_fig(fig_calib, 'calibration', 7.4, 4.4)

# simulation recovery
load(paste0(path_out, 'snap_recovery.rda'))
rec = snap_recovery$snap_prob_recovered
truth = snap_recovery$delta_true
n_sim = nrow(rec)
t_sim = ncol(rec)
heat_df = data.frame(
	actor = rep(1:n_sim, times = t_sim),
	year = rep(1:t_sim, each = n_sim),
	prob = as.numeric(rec),
	planted = as.numeric(truth) == 1
)
heat_df = heat_df[heat_df$year > 1, ]
plant_df = heat_df[heat_df$planted, ]
fig_recov = ggplot(heat_df, aes(x = year, y = actor, fill = prob)) +
	geom_tile(color = 'grey90') +
	geom_point(data = plant_df, aes(x = year, y = actor), shape = 4,
		size = 2.4, stroke = 0.9, inherit.aes = FALSE) +
	scale_fill_gradient(low = 'white', high = 'black',
		name = 'Posterior Snap Probability', limits = c(0, 1)) +
	scale_y_continuous(breaks = NULL) +
	labs(x = 'Simulated Year', y = 'Simulated State') +
	theme_paper
save_fig(fig_recov, 'sim_recovery', 7, 5)

# simulation specificity
spec_df = rbind(
	data.frame(
		prob = as.numeric(snap_recovery$snap_prob_drift_only[, -1]),
		source = 'Drift-Only Data (No Snaps Planted)'
	),
	data.frame(
		prob = snap_recovery$pos_probs,
		source = 'At Planted Snaps'
	)
)
spec_df$source = factor(spec_df$source,
	levels = c('Drift-Only Data (No Snaps Planted)', 'At Planted Snaps'))
fig_spec = ggplot(spec_df, aes(x = prob)) +
	geom_histogram(bins = 30, boundary = 0) +
	facet_wrap(~ source, ncol = 1, scales = 'free_y') +
	labs(x = 'Posterior Snap Probability', y = 'Count') +
	theme_paper
save_fig(fig_spec, 'sim_specificity', 7, 4.5)

# domain comparison
dom_series = list(
	'Diplomatic (UN Voting)' = rupture_index$un_voting$index_by_year,
	'Economic (Trade)' = rupture_index$trade$index_by_year,
	'Coercion (Sanctions: Sender)' = rupture_index$sanctions$index_by_year,
	'Coercion (Sanctions: Receiver)' = rupture_index$sanctions$receiver_index_by_year,
	'Coercion (Sanctions: Sender or Receiver)' = rupture_index$sanctions$combined_index_by_year
)
dom = do.call(rbind, lapply(names(dom_series), function(rel){
	idx = dom_series[[rel]]
	data.frame(year = as.integer(names(idx)), index = as.numeric(idx), relation = rel)
}))
dom = dom[is.finite(dom$index), ]
dom$relation = factor(dom$relation, levels = names(dom_series))
dom_2025 = dom[dom$year == 2025, ]
fig_dom = ggplot(dom[dom$year >= 2010, ], aes(x = year, y = index)) +
	geom_vline(xintercept = 2025, linetype = 'dashed', color = 'grey65', linewidth = 0.3) +
	geom_line() +
	geom_point(size = 0.8) +
	geom_point(data = dom_2025, color = '#d62728', size = 1.7) +
	facet_wrap(~ relation, ncol = 1) +
	scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
	labs(x = 'Year', y = 'Focal-Actor Rupture Index') +
	theme_paper
save_fig(fig_dom, 'domains', 7, 8)
