
# 
entropy = function(x, na.rm=T)  {x = x[!is.na(x)]; x =x[x!=0]; p=x/sum(x); if (length(x) != 0) {return(-sum(p*log(p)))} else {NA}}

#
nentropy = function(x, na.rm=T)  {y = x[!is.na(x)]; y =y[y!=0]; p=y/sum(y); if (length(y) != 0) {return(-sum(p*log(p))/log(length(x)))} else {NA}}

#
tau = function(x, na.rm=T) {return(sum(1-x/max(x))/(length(x)-1))}

# 
cv = function(x, na.rm=T) {sd(x, na.rm=na.rm)/mean(x, na.rm=na.rm)}

# 
logit = function(x, base=exp(1)) {return(log(x/(1-x), base=base))}

# range
range = function(x, na.rm=T) {return(max(x, na.rm=na.rm) - min(x, na.rm=na.rm))}

# Dynamic range
dynRange = function(x, na.rm=T) {
	#if (is.character(x) || is.factor(x)) {cat('char\n');return(NA)}
	x = as.numeric(x);
	x = replace(x, x==0, NA);
	if (sum(!is.na(x)) < 2) {
		return(NA)
	} else{
		return(log10(max(x, na.rm=T)) - log10(min(x, na.rm=T)))
	}
}


# Shuffle data.frame
shuffleData = function(df) {
	s = matrix(sample(unlist(m), ncol(m)*nrow(m)), ncol=ncol(m), nrow=nrow(m))
	rownames(s) <- rownames(df)
	colnames(s) <- colnames(df)
	return(as.data.frame(s))
}

# Extract the legend from ggplot
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}



#~~~~~~~~~~~~~~~~~~~~~~
# ggplot2 custom theme
#~~~~~~~~~~~~~~~~~~~~~~

# Modified from theme_bw() from this site http://sape.inf.usi.ch/quick-reference/ggplot2/themes

#theme_minimal <- function(base_size = 12) {
#  library(grid)
#  structure(list(
#    axis.line =         element_blank(),
#    axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, vjust = 1),
#    axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, hjust = 1),
#    axis.ticks =        element_line(colour = "black", size = 0.2),
#    axis.title.x =      element_text(size = base_size, vjust = 1),
#    axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
#    axis.ticks.length = unit(0.3, "lines"),
#    axis.ticks.margin = unit(0.5, "lines"),
# 
#    legend.background = element_rect(colour=NA), 
#    legend.key =        element_rect(colour = "grey80"),
#    legend.key.size =   unit(1.2, "lines"),
#    legend.text =       element_text(size = base_size * 0.8),
#    legend.title =      element_text(size = base_size * 0.8, face = "bold", hjust = 0),
#    legend.position =   "right",
# 
#    panel.background =  element_rect(fill = "white", colour = NA), 
#    panel.border =      element_rect(fill = NA, colour="grey50"), 
#    panel.grid.major =  element_line(colour = "grey90", size = 0.2),
#    panel.grid.minor =  element_line(colour = "grey98", size = 0.5),
#    panel.margin =      unit(0.25, "lines"),
# 
#    strip.background =  element_rect(fill = "grey80", colour = "grey50"), 
#    strip.text.x =      element_text(size = base_size * 0.8),
#    strip.text.y =      element_text(size = base_size * 0.8, angle = -90),
# 
#    plot.background =   element_rect(colour = NA),
#    plot.title =        element_text(size = base_size * 1.2),
#    plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")
#  ), class = "options")
#}


