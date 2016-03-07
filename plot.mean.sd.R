#!/usr/bin/env Rscript

options(stringsAsFactors=F)

# ==================
# DEBUG OPTIONS
# ==================

opt = list()
opt$output = "detected.features.pdf"
opt$merge_mdata_on = "labExpId"
opt$base_size = 16

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input_matrix"), default="stdin", 
    help="the matrix you want to analyze [default=%default]"),

make_option(c("-o", "--output"), default="detected.features.pdf",
    help="output file name. [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="input file has header [default=%default]"),

make_option(c('--log10'), action="store_true", default=FALSE,
	help="Apply the log to the data [default=%default]"),

make_option(c("-x", "--x_index"), default=1, type="character",
	help="column index (or indeces) for the x axis. [default=%default]"),

make_option(c("-y", "--y_index"), default=2, type="integer",
	help="column index for the y axis. [default=%default]"),

make_option(c("-F", "--fun"), type="character", default="mean_sdl",
	help="function to aggregate [default=%default]"),

make_option(c("-G", "--geom"), type="character", default="pointrange",
	help="function to aggregate: pointrange | bar | crossbar [default=%default]"),

make_option(c("-c", "--color"), type="character", default="orange",
	help="color of the point or bar [default=%default]"),

make_option(c("--fill"), type="character", default="orange",
	help="fill of the point or bar. [default=%default]"),
#	help="fill of the point or bar. If an integer is provided, is the index of the column for the fill factor [default=%default]"),

make_option(c("-f", "--facet_by"),  
	help="column index to facet by"),

make_option(c("--facet_scale"), default="fixed",
	help="facet scale: fixed | free | free_x | free_y"),

make_option(c("--y_title"), default="Percentage of detected features",
	help="title for the y axis [default=%default]"),

make_option(c("-B", "--base_size"), default=16,
	help="font base size [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help="width of the plot in inches [default=%default]"),

make_option(c("-H", "--height"), default=5,
	help="height of the plot in inches [default=%default]")

)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

##------------
## LIBRARIES
##------------
 
cat("Loading libraries... ")
#suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
source("~/R/functions.R")
cat("DONE\n\n")


#########
# BEGIN #
#########

# read input matrix (decide if including the header or not)
if (opt$input_matrix == "stdin") {
	m = read.table(file("stdin"), h=opt$header, sep="\t")
} else {
	m = read.table(opt$input_matrix, h=opt$header, sep="\t")
}

df = m


# Define x and y according to user indeces
#x = colnames(df)[opt$x_index]
x_indeces = as.numeric(strsplit(opt$x_index, split=",")[[1]])
df$x_labels = apply(df[x_indeces], 1, paste, collapse=',')
x = 'x_labels'
y = colnames(df)[opt$y_index]

# Apply the log if needed
if (opt$log10) {df[y] <- log10(df[y])}

# Apply the function outside ggplot
#tmp = as.data.frame(as.list(aggregate(formula, df, smedian.hilow)))
#tmp = rename(tmp, replace=c("V4.Median" = "y", "V4.Lower"="ymin", "V4.Upper"="ymax"))


# sort by mean
formula = as.formula(sprintf("%s~%s", y, x))
lev = aggregate(formula, df, mean, na.rm=T)
lev = lev[order(lev[,y], decreasing=T),x]
df[x] = factor(df[,x], levels=lev)

# ~~~~~~~~~
# PLOT
# ~~~~~~~~~

theme_set(theme_bw(base_size=opt$base_size))

gp = ggplot(df, aes_string(x=x, y=y))
if (opt$geom == 'crossbar') {
	gp = gp + stat_summary(fun.data=opt$fun, mult=1, shape=15, size=1, width=0.4, color=opt$color, fill=opt$fill, geom=opt$geom)
} else {
	gp = gp + stat_summary(fun.data=opt$fun, mult=1, shape=15, size=1, width=0.4, color=opt$color, geom="errorbar")
	gp = gp + stat_summary(fun.data=opt$fun, mult=1, shape=15, size=1, color=opt$color, fill=opt$fill, geom=opt$geom)
}
gp = gp + theme(axis.text.x=element_text(angle=45, hjust=1))
gp = gp + labs(y=opt$y_title, x="")


if (!is.null(opt$facet_by)) {
	facet = as.formula(sprintf("~%s", colnames(df)[as.numeric(opt$facet_by)]))
	gp = gp + facet_wrap(facet, nrow=1, scale=opt$facet_scale)
}


ggsave(opt$output, h=opt$height, w=opt$width)

q(save="no")


