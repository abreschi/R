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

make_option(c("-b", "--binwidth"), type="double", 
	help="Specify binwidth. Leave empty for default")

)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
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


avg = mean(df[,1], na.rm=TRUE)
med = median(df[,1], na.rm=TRUE)


# GGPLOT

theme_set(theme_bw(base_size=20))

gp = ggplot(df, aes_string(x = colnames(df)[1]))
if (is.null(opt$binwidth)) {
	gp = gp + geom_histogram(fill=opt$fill, color=opt$color, right=TRUE, include.lowest=TRUE)
} else {
	gp = gp + geom_histogram(fill=opt$fill, color=opt$color, right=TRUE, include.lowest=TRUE, binwidth=opt$binwidth)
}
gp = gp + geom_point(aes(x=avg, y=0), size=2)
gp = gp + geom_vline(xintercept=med, linetype=2)
gp = gp + labs(y=opt$y_title, x=opt$x_title)

ggsave(opt$output, h=5, w=5)

# EXIT
quit(save='no')
