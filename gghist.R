#!/usr/bin/env Rscript


options(stringsAsFactors=FALSE)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin [default=%default]"),

make_option(c("-o", "--output"), default="gghist.out.pdf",
	help="Output file name [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

make_option(c("--scale_x_log10"), action="store_true", default=FALSE,
	help="log10-transform x scale [default=%default]"),

make_option(c("--scale_y_log10"), action="store_true", default=FALSE,
	help="log10-transform y scale [default=%default]"),

make_option(c("--y_title"), type="character", default="count",
	help="Title for the y axis [default=%default]"),

make_option(c("--x_title"), type="character", default="",
	help="Title for the x axis [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]"),

make_option(c("-f", "--fill"), default="aquamarine",
	help="choose the color which you want to fill the histogram with"),

make_option(c("-c", "--color"), default="grey",
	help="choose the color which you want to contour the histogram with"),

make_option(c("-F", "--fill_by"), type='numeric',
	help="the column index with the factor to fill by. Leave empty for no factor."),

make_option(c("-W", "--width"), default=7,
	help="width of the plot in inches. [default=%default]"),

make_option(c("-b", "--binwidth"), type="double", 
	help="Specify binwidth. Leave empty for default")

)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list,
	description = "Reads the values on the first column and outputs a histogram"
)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

#------------
# LIBRARIES
#------------ 

if (opt$verbose) {cat("Loading libraries... ")}
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


# ======== #
# BEGIN    #
# ======== #


# Read data

if (opt$input == "stdin") {
	m = read.table(file("stdin"), h=opt$header) 
} else {
	m = read.table(opt$input, h=opt$header)
}

df = m


if (!is.null(opt$fill_by)) {
	fill_by = colnames(df)[opt$fill_by]
}


# GGPLOT

theme_set(theme_bw(base_size=20))

gp = ggplot(df, aes_string(x = colnames(df)[1]))
if (is.null(opt$binwidth)) {
	gp = gp + geom_histogram(fill=opt$fill, color=opt$color, right=TRUE, include.lowest=TRUE)
} else {
	gp = gp + geom_histogram(fill=opt$fill, color=opt$color, right=TRUE, include.lowest=TRUE, binwidth=opt$binwidth)
}

if (!is.null(opt$fill_by)) {
	gp = gp + geom_histogram(aes_string(fill=fill_by))
}

gp = gp + theme(axis.text.x=element_text(angle=45, hjust=1))

if (!is.character(df[,1])) {
	avg = mean(df[,1], na.rm=TRUE)
	med = median(df[,1], na.rm=TRUE)
	gp = gp + geom_point(aes(x=avg, y=0), size=2)
	gp = gp + geom_vline(xintercept=med, linetype=2)
	gp = gp + theme(axis.text.x=element_text(angle=0, hjust=0.5))
}


if (opt$scale_x_log10) {
	gp = gp + scale_x_log10()
}

if (opt$scale_y_log10) {
	gp = gp + scale_y_log10()
}

gp = gp + labs(y=opt$y_title, x=opt$x_title)

ggsave(opt$output, h=5, w=opt$width)

# EXIT
quit(save='no')
