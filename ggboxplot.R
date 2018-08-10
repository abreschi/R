#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
set.seed(123)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), type="character", default='stdin',
	help="tab-separated file. Can be stdin [default=%default]"),

make_option(c("-o", "--output"), default="ggboxplot.pdf", 
	help="output file name with extension [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE, 
	help="The file has header [default=%default]"),

make_option(c("--xy"), type='character', default="1,2",
	help="Column indeces of the x and y axes, respectively, comma-separated [default=%default]"),

make_option(c("-c", "--color_by"), type='numeric',
	help="Column index with the color factor. Leave empty for no color"),

make_option(c("-f", "--fill_by"), type='numeric', 
	help="Column index of the fill factor, Leave empty for no fill"),

make_option(c("-r", "--layers"), default="boxplot",
	help="Choose layers to add on the plot (comma-separated): < boxplot | violin | jitter > [default=%default]"),

make_option(c("--facet_by"), type="integer", 
	help="column index to facet"),

make_option(c("--facet_scale"), default="fixed",
	help="Scale for faceting [default=%default]"),

make_option(c("--facet_nrow"), type="integer",
	help="Number of rows in facet wrap [default=%default]"),

make_option(c("--log"), action="store_true", default=FALSE,
	help="apply the log to the y-axis [default=%default]"),

make_option(c("-p", "--pseudocount"), type="numeric", default=1e-03,
	help="Pseudocount for the log on the y-axis [default=%default]"),

make_option(c("-t", "--title"), type="character", 
	help="Main title for the plot. Leave emtpy for no title"),

make_option(c("--y_title"),  
	help="title for the y-axis. Leave empty for no title"),

make_option(c("--x_title"), 
	help="title for the x-axis. Leave empty for no title"),

make_option(c("--y_limits"), 
	help="Specify limits for the y axis, e.g. \"\\-1,1\". Escape character for negative numbers [default=%default]"),

make_option(c("-P", "--palette"), 
	help='File with colors for the lines. Leave empty to use even color spacing'),

make_option(c("-a", "--annotate"), action='store_true', default=FALSE,
	help="Print median over boxplots [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help='Width of the plot in inches [default=%default]'),

make_option(c("-H", "--height"), default=5,
	help='Height of the plot in inches [default=%default]'),

make_option(c("--x_order"), default=NULL,
	help="File with the order of x labels [default=%default]"),

make_option(c("-v", "--verbose"), action='store_true', default=FALSE,
	help="Verbose output [default=%default]")

#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list = option_list,
	description = "From a column file, plot a boxplot or a violin plot"
)


arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


##------------
## LIBRARIES
##------------
 
if (opt$verbose) {cat("Loading libraries... ")}
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


# ========
# BEGIN
# ========


# Read data
if (opt$input == "stdin") {input=file('stdin')} else {input=opt$input}
m = read.table(input, h=opt$header, sep="\t")
if(opt$verbose) {print(head(m))}

# Read the axes
xy = strsplit(opt$xy, ",")[[1]]
x_col = as.numeric(xy[1])
y_col = as.numeric(xy[2])
x = colnames(m)[x_col]
y = colnames(m)[y_col]
if (!is.null(opt$color_by)) {color_by = colnames(m)[opt$color_by]} else {color_by=NULL}
if (!is.null(opt$fill_by)) {fill_by = colnames(m)[opt$fill_by]} else {fill_by=NULL}
if (!is.null(opt$facet_by)) {facet_col = colnames(m)[opt$facet_by]}

# Convert to log if asked
if (opt$log) {
	m[,y] <- log10(m[,y]+opt$pseudocount)
}

# Read palette
if (!is.null(opt$palette)) {
	palette = read.table(opt$palette, h=FALSE, comment.char='%', sep="\t")$V1
}


# # Significance
# levs = levels(as.factor(m[,x]))
# pv = sapply(1:(length(levs)-1), function(i) {
# 	distr1 = m[which(m[,x] == levs[i]), y];
# 	distr2 = m[which(m[,x] == levs[i+1]), y];
# 	return(suppressWarnings(ks.test(distr1, distr2)$p.value))
# 	}
# )
# sign_df = data.frame(x=levs[-length(levs)], text=format(pv, digits=2))

# Replace newlines if present in character columns

charCols = sapply(m, class) == "character"
for (i in colnames(m)[charCols]) {
	m[,i] <- gsub("\\\\n", "\n", m[,i])
}


# Compute the quantiles for the data based on grouping
if (opt$annotate) {
	medianFormula = as.formula(sprintf("%s~%s", y, x)) 
	if (!is.null(color_by) || !is.null(fill_by) || !is.null(facet_by)) {
		medianFormula = as.formula(sprintf("%s~%s", y, paste(x, color_by, fill_by, facet_col, sep="+")))
	}
	medians = aggregate(medianFormula, data=m, median)
	colnames(medians)[length(colnames(medians))] = "median"
	medians$median <- round(medians$median, digits=2)
	uppQuartiles = aggregate(medianFormula, data=m, quantile, 0.75)
	colnames(uppQuartiles)[length(colnames(uppQuartiles))] = "uppQuartile"
	text = merge(medians, uppQuartiles, by=intersect(colnames(medians), colnames(uppQuartiles)))
}



#~~~~~~~~~~~~
# GGPLOT
#~~~~~~~~~~~~

theme_set(theme_bw(base_size=20))
theme_update(
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	panel.spacing = unit(0.2, "lines"),
	legend.key = element_blank(),
	plot.title = element_text(vjust=1),
	axis.title.y = element_text(vjust=1.4,angle=90)
)


mapping = aes_string(color=color_by, fill=fill_by)

alpha = 0.9

# Read x order
if (!is.null(opt$x_order)) {
	x_order = read.table(opt$x_order, h=F)[,1]
	m[,x] = factor(m[,x], levels=x_order)
}

# plot 
gp = ggplot(m, aes_string(x=x, y=y))

layers = strsplit(opt$layers, split=",")[[1]]
for (l in layers) {
	geom_params = list()

	geom_params$alpha = alpha
	geom_params$size = 1

	if (l == "boxplot") {
		geom = "boxplot"
		stat = "boxplot"
		position = "dodge"
		if ("jitter" %in% layers) {
			geom_params$outlier.shape = NA
			mapping = aes_string()
		}
	}
	if (l == "jitter") {
		geom = "point"
		stat = "identity"
		position = "jitter"
		mapping = aes_string(color=color_by, fill=fill_by)
	}
	plotLayer <- layer(
		geom = geom,
		stat = stat,
		params = geom_params,
		mapping = mapping,
		position = position
	)
	gp = gp + plotLayer
}

if (opt$annotate) { ##### TO DO: write it as layers!!! 
	gp = gp + geom_text(data=text, aes_string(x=x, y="uppQuartile", label="median"), vjust=-2)
}

if (!is.null(opt$y_title)) {
	opt$y_title = gsub("\\\\n", "\n", opt$y_title)
}

if (!is.null(opt$x_title)) {
	opt$x_title = gsub("\\\\n", "\n", opt$x_title)
}

if (!is.null(opt$title)) {
	opt$title = gsub("\\\\n", "\n", opt$title)
}

gp = gp + labs(title=opt$title, y=opt$y_title, x=opt$x_title)
gp = gp + theme(axis.text.x=element_text(angle=45, hjust=1))

# Fill scale
if (!is.null(opt$palette)) {
	gp = gp + scale_fill_manual(values=palette)
	gp = gp + guides(fill=guide_legend(ncol=1))
} else {
	gp = gp + scale_fill_hue()
}
# Color scale
if (!is.null(opt$palette)) {
	gp = gp + scale_color_manual(values=palette)
	gp = gp + guides(color=guide_legend(ncol=1))
} else {
	gp = gp + scale_color_hue()
}

y_limits = NULL	
if (!is.null(opt$y_limits)) {
	opt$y_limits = gsub("\\", "", opt$y_limits, fixed=TRUE)
	y_limits = as.numeric(strsplit(opt$y_limits, ",")[[1]])
} 
gp = gp + scale_y_continuous(limits=y_limits)

# Read facet
if (!is.null(opt$facet_by)) {
	facet_col = colnames(m)[opt$facet_by]
	facet_form = as.formula(sprintf("~%s", facet_col))
	gp = gp + facet_wrap(facet_form, scales=opt$facet_scale, nrow=opt$facet_nrow)
}

gp = gp + guides(fill=guide_legend(ncol=1))
gp = gp + guides(color=guide_legend(ncol=1))

# adjust size for long x labels 
gt = ggplot_gtable(ggplot_build(gp))
gt$widths[1] = unit(as.numeric(convertUnit(unit(1, 'strwidth', 
	levels(as.factor(m[[x]]))[1] ), 'in')) / sqrt(2),'in') - gt$widths[2]

pdf(opt$output, h=opt$height, w=opt$width)
grid.draw(gt)
dev.off()
#ggsave(opt$output, h=opt$height, w=opt$width)

# EXIT
quit(save='no')
