#!/usr/bin/env Rscript

##------------
## LIBRARIES
##------------
 
cat("Loading libraries... ")
#suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
source("~/R/functions.R")
cat("DONE\n\n")


options(stringsAsFactors=F)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


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

option_list <- list(

make_option(c("-i", "--input_matrix"), default="stdin", 
    help="the matrix you want to analyze [default=%default]"),

make_option(c("-o", "--output"), default="detected.features.pdf",
    help="output file name. [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="input file has header"),

make_option(c("-x", "--x_index"), default=1, type="integer",
	help="column index for the x. [default=%deafult]"),

make_option(c("-y", "--y_index"), default=2, type="integer",
	help="column index for the y. [default=%deafult]"),

#make_option(c("-m", "--metadata"), 
#	help="metadata index file"),
#
#make_option(c("--merge_mdata_on"), default="labExpId",
#	help="field to merge the metadata on [default=%default]"),

make_option(c("-f", "--facet_by"), type="integer", 
	help="column index to facet by"),

make_option(c("--mean_by"), default="labExpId",
	help="metadata field to average by. [default=%default]"),

make_option(c("--y_title"), default="Percentage of detected features",
	help="title for the y axis [default=%default]"),

make_option(c("-B", "--base_size"), default=16,
	help="font base size [default=%default]")

)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#########
# BEGIN #
#########

# read input matrix (decide if including the header or not)
if (opt$input_matrix == "stdin") {
	m = read.table(file("stdin"), h=opt$header)
} else {
	m = read.table(opt$input_matrix, h=opt$header)
}

# read metadata, if provided
if (!is.null(opt$metadata)) {
	mdata = read.table(opt$metadata, h=T, sep="\t")
	mdata[,opt$merge_mdata_on] <- gsub(",", ".", mdata[,opt$merge_mdata_on])
	mdata = unique(mdata[c(opt$facet_by, opt$merge_mdata_on, opt$mean_by)])
	df = merge(m, mdata, by.x="V1", by.y=opt$merge_mdata_on)
}

df = m

# plot
theme_set(theme_bw(base_size=opt$base_size))

x = colnames(df)[opt$x_index]
y = colnames(df)[opt$y_index]

gp = ggplot(df, aes_string(x=x, y=y))
gp = gp + stat_summary(fun.data="mean_sdl", mult=1, shape=15, size=1, color="blue")
gp = gp + theme(axis.text.x=element_text(angle=45, hjust=1))
gp = gp + labs(y=opt$y_title, x="")

if (!is.null(opt$facet_by)) {
	facet = as.formula(sprintf("~%s", colnames(df)[opt$facet_by]))
	gp = gp + facet_wrap(facet, nrow=1)
}


ggsave(opt$output, h=5, w=7)

q(save="no")


