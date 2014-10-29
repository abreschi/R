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

#make_option(c("-f", "--facet"), type="integer", help="column index to facet"),

make_option(c("-t", "--title"), type="character", 
	help="Main title for the plot. Leave emtpy for no title"),

make_option(c("--y_title"), default="norm_read_density", 
	help="title for the y-axis [default=%default]"),

make_option(c("--x_title"), default="position", 
	help="title for the x-axis [default=%default]"),

make_option(c("-P", "--palette"), 
	help='File with colors for the lines. Leave empty to use even color spacing'),

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
m = read.table(input, h=opt$header)
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
	panel.grid.major = element_blank()
)

x = colnames(m)[opt$x_col]
y = colnames(m)[opt$y_col]
if (!is.null(opt$color_by)) {color_by = colnames(m)[opt$color_by]} else {color_by=NULL}
if (!is.null(opt$linetype_by)) {linetype_by = colnames(m)[opt$linetype_by]} else {linetype_by=NULL}
alpha = 1

# plot 
gp = ggplot(m, aes_string(x=x, y=y)) 
gp = gp + geom_line(aes_string(color=color_by, linetype=linetype_by), alpha=alpha, size=1)
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


ggsave(opt$output, h=5, w=7)

# EXIT
quit(save='no')
