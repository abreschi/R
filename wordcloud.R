#!/usr/bin/env Rscript


options(stringsAsFactors=FALSE)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

make_option(c("-o", "--output"), default="wordcloud.pdf",
	help="Output file name [default=%default]"),

make_option(c("-w", "--words"), default=1,
	help="Index of the column with the words [default=%default]"),

make_option(c("-e", "--weights"), default=2, type="integer",
	help="Index of the column with weights [default=%default]"),

make_option(c("-C", "--color_by"), type="integer", 
	help="Index of the factor for coloring"),

make_option(c("-P", "--palette"), default="/users/rg/abreschi/R/palettes/Spectral.11.txt",
        help='File with colors [default=%default]'),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")

#make_option(c("--position"), default='dodge',
#	help="Position for histogram [default=%default]"),
#
#make_option(c("--scale_x_log10"), action="store_true", default=FALSE,
#	help="log10-transform x scale [default=%default]"),
#
#make_option(c("--scale_y_log10"), action="store_true", default=FALSE,
#	help="log10-transform y scale [default=%default]"),
#
#make_option(c("--y_title"), type="character", default="count",
#	help="Title for the y axis [default=%default]"),
#
#make_option(c("--x_title"), type="character", default="",
#	help="Title for the x axis [default=%default]"),
#
#make_option(c("-f", "--fill"), default="aquamarine",
#	help="choose the color which you want to fill the histogram with"),
#
#make_option(c("-c", "--color"), default="white",
#	help="choose the color which you want to contour the histogram with"),
#
#make_option(c("-F", "--fill_by"), type='numeric',
#	help="the column index with the factor to fill by. Leave empty for no factor."),
#
#make_option(c("--facet_by"), type='numeric',
#	help="the column index with the factor to facet by. Leave empty for no factor."),
#
#make_option(c("-W", "--width"), default=7,
#	help="width of the plot in inches. [default=%default]"),
#
#make_option(c("-b", "--binwidth"), type="double", 
#	help="Specify binwidth. Leave empty for default")

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
#suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(wordcloud))
if (opt$verbose) {cat("DONE\n\n")}


# ======== #
# BEGIN    #
# ======== #

set.seed(123)

# Read data
if (opt$input == "stdin") {input=file("stdin")} else {input=opt$input}
m = read.table(input, h=opt$header, sep="\t") 

df = m

words = df[,opt$words] 
freq = df[,opt$weights]

if (!is.null(opt$color_by)) {
	palette = as.character(read.table(opt$palette, comment.char="%", h=F)$V1)
	color_by = df[,opt$color_by]
	colors = palette[as.factor(color_by)]
}


pdf(opt$output, width=7, height=7)
wordcloud(
	words = words, 
	freq = freq, 
	min.freq = -1,
	colors = colors,
	ordered.colors = TRUE
)
dev.off()

q(save='no')
