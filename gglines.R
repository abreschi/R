#!/usr/bin/env Rscript

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), type="character", default='stdin',
	help="list of files with the bigWig profiles. Can be stdin [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE, 
	help="The file has header [default=%default]"),

make_option(c("-c", "--color_by"), type='numeric',
	help="Column index with the color factor. Leave empty for no color"),

make_option(c("-x", "--x_col"), type='numeric', default=1,
	help="Column index of the x axis [default=%default]"),

make_option(c("-y", "--y_col"), type='numeric', default=2,
	help="Column index of the y axis [default=%default]"),

#make_option(c("-L", "--labels"), help="list of labels with the labels of each file, commma-separated.\n
#They must be in the same order as the file list", type="character"),
make_option(c("-v", "--vertical_lines"), type='character',
	help="specify where you want the vertical lines [default=%default]"),

#make_option(c("-f", "--facet"), type="integer", help="column index to facet"),
make_option(c("-t", "--title"), type="character", 
	help="Main title for the plot. Leave emtpy for no title"),

make_option(c("-o", "--output"), help="output file name with extension [default=%default]", default="profile.pdf"),
make_option(c("--y_title"), help="title for the y-axis [default=%default]", default="norm_read_density"),
make_option(c("-P", "--palette"), 
	help='File with colors for the lines. Leave empty to use even color spacing'),

make_option(c("-v", "--verbose"), action='store_true', default=FALSE,
	help="Verbose output [default=%default]")

#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
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
#if (!is.null(opt$labels)) {labels = strsplit(opt$labels, ",")[[1]]}

if (!is.null(opt$vertical_lines)) {vertical_lines = as.numeric(strsplit(opt$vertical_lines, ",")[[1]])}


# Read data
if (opt$input == "stdin") {input=file('stdin')} else {input=opt$input}
m = read.table(input, h=opt$header)

# Read palette
if (!is.null(opt$palette)) {
	palette = read.table(opt$palette, h=FALSE, comment.char='%')$V1
}


#~~~~~~~~~~~~
# GGPLOT
#~~~~~~~~~~~~

theme_set(theme_bw(base_size=20))
theme_update(panel.grid.minor=element_blank())

x = colnames(m)[opt$x_col]
y = colnames(m)[opt$y_col]
color_by = colnames(m)[opt$color_by]
alpha = 1

gp = ggplot(m, aes_string(x=x, y=y)) 
#gp = gp + geom_line(aes_string(group=color_by, color=color_by), alpha=0.5, size=1)
gp = gp + geom_line(aes_string(color=color_by), alpha=alpha, size=1)
# Color scale
if (!is.null(opt$palette)) {
	gp = gp + scale_color_manual(values=palette)
} else {
	gp = gp + scale_color_hue()
}
#gp = gp + geom_vline(xintercept=vertical_lines, color="grey", linetype="longdash")
#gp = gp + geom_vline(xintercept=NULL)
gp = gp + labs(title=opt$title, y=opt$y_title)

ggsave(opt$output, h=5, w=7)

# EXIT
quit(save='no')
