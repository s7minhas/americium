# prep ame circle plots

load_pkg('ggrepel')

get_data_for_circ = function(
	y,
	u = NULL,
	v = NULL,
	row_names = rownames(y),
	col_names = colnames(y),
	vscale = 0.8,
	remove_isolates = TRUE,
	u_label = 'U',
	v_label = 'V'
){
	v_missing = is.null(v)
	u_missing = is.null(u)
	
	if(u_missing){
		a = rowMeans(y, na.rm = TRUE)
		b = colMeans(y, na.rm = TRUE)
		y0 = y
		y0[is.na(y)] = outer(a, b, '+')[is.na(y)]
		y0 = y0 - mean(y0)
		
		if(!all(y == t(y), na.rm = TRUE)){
			y_svd = svd(y0)
			u = y_svd$u[, 1:2]
			v = y_svd$v[, 1:2]
			u_norm = sqrt(rowSums(u^2))
			v_norm = sqrt(rowSums(v^2))
			u = diag(1 / u_norm) %*% u
			v = diag(1 / v_norm) %*% v * vscale
		}
		
		if(all(y == t(y), na.rm = TRUE)){
			y_eigen = eigen(y0)
			keep = which(
				abs(y_eigen$values) >= sort(abs(y_eigen$values), decreasing = TRUE)[2]
			)[1:2]
			u = y_eigen$vectors[, keep]
			u_norm = sqrt(rowSums(u^2))
			u = diag(1 / u_norm) %*% u
			v_norm = u_norm
			v = u
		}
	}
	
	if(!u_missing){
		if(v_missing){
			v = u
			vscale = 1
		}
		u_norm = sqrt(rowSums(u^2))
		v_norm = sqrt(rowSums(v^2))
		u = diag(1 / u_norm) %*% u
		v = diag(1 / v_norm) %*% v * vscale
	}
	
	rownames(y) = row_names
	colnames(y) = col_names
	row_activity = rowSums(abs(y), na.rm = TRUE)
	col_activity = colSums(abs(y), na.rm = TRUE)
	links = which(y != 0, arr.ind = TRUE)
	
	u_data = data.frame(u * 1.2)
	u_data$actor = rownames(y)
	u_data$point_size = 0
	u_data$point_size[row_activity > 0] = u_norm[row_activity > 0]^3
	if(remove_isolates){
		u_data = u_data[u_data$point_size > 0, ]
	}
	
	if(!v_missing){
		v_data = data.frame(v * 1.2)
		v_data$actor = colnames(y)
		v_data$point_size = 0
		v_data$point_size[col_activity > 0] = v_norm[col_activity > 0]^3
		if(remove_isolates){
			v_data = v_data[v_data$point_size > 0, ]
		}
		
		u_data$effect = u_label
		v_data$effect = v_label
		u_data = rbind(u_data, v_data)
		u_data$effect = factor(u_data$effect, levels = c(u_label, v_label))
	}
	
	list(
		node_data = u_data,
		u = u,
		v = v,
		links = links
	)
}

# plot ame circles
gg_circ = function(
	y,
	u = NULL,
	v = NULL,
	row_names = rownames(y),
	col_names = colnames(y),
	vscale = 0.8,
	point_range = c(2, 5),
	link_color = 'gray85',
	link_type = 'dotted',
	link_size = 0.5,
	force = 1,
	max_iter = 3000,
	remove_isolates = TRUE,
	u_label = 'U',
	v_label = 'V',
	show_links = TRUE,
	show_labels = TRUE,
	show_text = FALSE,
	show_points = TRUE,
	...
){
	has_v = !is.null(v)
	circ_data = get_data_for_circ(
		y = y,
		u = u,
		v = v,
		row_names = row_names,
		col_names = col_names,
		vscale = vscale,
		remove_isolates = remove_isolates,
		u_label = u_label,
		v_label = v_label
	)
	node_data = circ_data$node_data
	u = circ_data$u
	v = circ_data$v
	links = circ_data$links
	
	if(has_v){
		plot_data = ggplot(node_data, aes(x = X1, y = X2, color = effect))
	}
	if(!has_v){
		plot_data = ggplot(node_data, aes(x = X1, y = X2))
	}
	
	if(show_links && nrow(links) > 0){
		for(i in seq_len(nrow(links))){
			plot_data = plot_data +
				geom_segment(
					x = u[links[i, 1], 1] * 1.2,
					y = u[links[i, 1], 2] * 1.2,
					xend = v[links[i, 2], 1] * 1.2,
					yend = v[links[i, 2], 2] * 1.2,
					color = link_color,
					linetype = link_type,
					size = link_size
				)
		}
	}
	
	if(show_points){
		plot_data = plot_data + geom_point()
	}
	if(show_labels){
		plot_data = plot_data +
			geom_label_repel(
				aes(label = actor, size = point_size, ...),
				force = force,
				max.iter = max_iter
			)
	}
	if(show_text){
		plot_data = plot_data +
			geom_text_repel(
				aes(label = actor, size = point_size, ...),
				force = force,
				max.iter = max_iter
			)
	}
	
	plot_data +
		scale_size(range = point_range) +
		theme_bw() +
		theme(
			legend.position = 'top',
			axis.ticks = element_blank(),
			axis.title = element_blank(),
			axis.text = element_blank(),
			panel.border = element_blank(),
			strip.background = element_rect(fill = 'black', color = 'black'),
			strip.text.x = element_text(color = 'white', hjust = 0),
			strip.text.y = element_text(color = 'white', hjust = 0)
		)
}
