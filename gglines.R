#!/usr/bin/env Rscript

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), type="character", default='stdin',
	help="tab-separated file. Can be stdin [default=%default]"),

make_option(c("-o", "--output"), default="profile.pdf", 
	help="output file name with extension [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE, 
	help="The file has header [default=%default]"),

make_option(c("--group"), type="numeric",
	help="Column index with the group factor [default=%default]"),

make_option(c("-c", "--color_by"), type='numeric',
	help="Column index with the color factor. Leave empty for no color"),

make_option(c("-L", "--linetype_by"), type='numeric', 
	help="Column index of the linetype factor, Leave empty for no linetype"),

make_option(c("-x", "--x_col"), type='numeric', default=1,
	help="Column index of the x axis [default=%default]"),

make_option(c("-y", "--y_col"), type='numeric', default=2,
	help="Column index of the y axis [default=%default]"),

make_option(c("-V", "--vertical_lines"), type='character', 
	help="specify where you want the vertical lines [default=%default]"),

make_option(c("-a", "--alpha"), type="numeric", default=1,
	help="Set transparency [default=%default]"),

#make_option(c("-f", "--facet"), type="integer", help="column index to facet"),

make_option(c("-t", "--title"), type="character", 
	help="Main title for the plot. Leave emtpy for no title"),

make_option(c("--y_title"), default="norm_read_density", 
	help="title for the y-axis [default=%default]"),

make_option(c("--x_title"), default="position", 
	help="title for the x-axis [default=%default]"),

make_option(c("--scale_x_log10"), action="store_true", default=FALSE,
	help="Change the x axis to log10 [default=%default]"),

make_option(c("--x_limits"), default=NULL,
	help="Limits for the x axis [default=%default]"),

make_option(c("-P", "--palette"), 
	help='File with colors for the lines. Leave empty to use even color spacing'),

make_option(c("-H", "--height"), default=5,
	help="Height of the plot in inches [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help="Width of the plot in inches [default=%default]"),

make_option(c("-v", "--verbose"), action='store_true', default=FALSE,
	help="Verbose output [default=%default]")

#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list = option_list,
	description = "From a column file, plot a column vs another as lines"
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
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


# ========
# BEGIN
# ========



# Read data
if (opt$input == "stdin") {input=file('stdin')} else {input=opt$input}
m = read.table(input, h=opt$header, sep="\t")
if(opt$verbose) {print(head(m))}

# Read palette
if (!is.null(opt$palette)) {
	palette = read.table(opt$palette, h=FALSE, comment.char='%')$V1
}

# Read coloring factor
if (!is.null(opt$color_by)) {
	if (ncol(m)<opt$color_by) {
		cat("ERROR: color factor index out of range\n")
		q(save='no')
	}
}


#~~~~~~~~~~~~
# GGPLOT
#~~~~~~~~~~~~

theme_set(theme_bw(base_size=20))
theme_update(
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	legend.key = element_blank()
)

x = colnames(m)[opt$x_col]
y = colnames(m)[opt$y_col]



alpha = opt$alpha

# -------------------------
# geom_params and mapping
# -------------------------

geom_params = list()
mapping = list()


# >> color <<

if (is.null(opt$color_by)) {
	geom_params$color = "black"
} else {
	mapping = modifyList(mapping, aes_string(color=colnames(m)[opt$color_by]))
}


# >> linetype <<

if (is.null(opt$linetype_by)) {
	geom_params$linetype = 1
} else {
	mapping = modifyList(mapping, aes_string(linetype=colnames(m)[opt$linetype_by]))
}


# >> group <<

if (!is.null(opt$group)) {
	group = colnames(m)[opt$group]
	mapping = modifyList(mapping, aes_string(group=colnames(m)[opt$group]))
	m[,group] <- gsub("\\\\n", "\n", m[,group])
}

geom_params$alpha = alpha
geom_params$size = 1

class(mapping) <- "uneval"

# -----------
# lineLayer
# -----------

lineLayer <- layer(
	geom = "line",
	params = geom_params,
#	geom_params = geom_params,
	mapping = mapping,
	stat = "identity",
	position = "identity"
)


# plot 
gp = ggplot(m, aes_string(x=x, y=y)) 

# Add line layer
gp = gp + lineLayer

# Add the labels
gp = gp + labs(title=opt$title, y=opt$y_title, x=opt$x_title)

# Color scale
if (!is.null(opt$palette)) {
	gp = gp + scale_color_manual(values=palette)
} else {
	gp = gp + scale_color_hue()
}

# Vertical lines
if (!is.null(opt$vertical_lines)) {
	vertical_lines = as.numeric(strsplit(opt$vertical_lines, ",")[[1]])
	gp = gp + geom_vline(xintercept=vertical_lines, color="grey", linetype="longdash")
}

# Scale x log10
if (opt$scale_x_log10) {
	gp = gp + scale_x_log10()
}

# Read x limits
if (!is.null(opt$x_limits)) {
	x_lim = as.numeric(strsplit(opt$x_limits, ",")[[1]])
	gp = gp + xlim(x_lim)
}

ggsave(opt$output, h=opt$height, w=opt$width)

# EXIT
quit(save='no')
