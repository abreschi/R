#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin [default=%default]"),

make_option(c("-o", "--output"), default="quantile.out.tsv",
	help="Output file name. <stdout> for printing on stdout [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

make_option(c("-m", "--method"), default="quantiles",
	help="Choose the way to bin the data: < quantiles | breaks > [default=%default]"),

make_option(c("-c", "--column"), default=1,
	help="The column with the values to divide in quantiles [default=%default]"),

make_option(c("-q", "--quantiles"), default=4, type="integer",
	help="Number of quantiles [default=%default]"),

make_option(c("-s", "--resolve_breaks"), default=FALSE, action="store_true",
	help="When breaks are not unique, don't crash but create fewer quantiles [default=%default]"),

make_option(c("-b", "--breaks"), default="-Inf,1,2,3,Inf", 
	help="Breaks for the intervals, comma-separated [default=%default]"), 

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")

)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

#------------
# LIBRARIES
#------------ 

if (opt$verbose) {cat("Loading libraries... ")}
#suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


# ============== #
# Functions      #
# ============== #

formatInterval = function(interval) {

	interval = as.character(interval)
	if (length(grep("-Inf", interval) != 0)) {
		outInterval = sprintf("<= %s", gsub("\\)", "", strsplit(interval, split=",")[[1]][2]))
		return(outInterval)
	}
	if (length(grep("Inf", interval) != 0)) {
		outInterval = sprintf("> %s", gsub("\\(", "", strsplit(interval, split=",")[[1]][1]))
		return(outInterval)
	}

#	left = gsub("[\\(\\]]", "", strsplit(interval, split=",")[[1]][1], perl=TRUE)
#	right = gsub("\\)", "", strsplit(interval, split=",")[[1]][2])
#	outInterval = sprintf("%s <= x < %s", left, right)
	outInterval = interval
	return(outInterval)
}


# ======== #
# BEGIN    #
# ======== #


# Read data

if (opt$input == "stdin") {
	m = read.table(file("stdin"), h=opt$header) 
} else {
	m = read.table(opt$input, h=opt$header)
}


# Quantiles

if (opt$method == "quantiles") {

quantile_header = sprintf("quantile_%s", colnames(m)[opt$column])
quant_index_header = sprintf("quant_index_%s", colnames(m)[opt$column])

if (opt$resolve_breaks) {
	breaks <- c(unique(quantile(m[,opt$column], probs=seq(0,1,1/opt$quantiles), na.rm=T)))
	m[,quantile_header] = cut(m[,opt$column], breaks, include.lowest=TRUE)
	m[,quant_index_header] = match(m[,quantile_header], levels(m[,quantile_header]))
} else {
	m[,quantile_header] = cut_number(m[,opt$column], opt$quantiles)
	m[,quant_index_header] = as.numeric(cut_number(m[,opt$column], opt$quantile))
}

}

# Intervals

if (opt$method == "breaks") {

interval_header = sprintf("interval_%s", colnames(m)[opt$column])
interval_index_header = sprintf("interval_index_%s", colnames(m)[opt$column])

breaks = as.numeric(strsplit(opt$breaks, split=",")[[1]])
m[,interval_header] = cut(m[,opt$column], breaks=breaks, include.lowest=TRUE, right=TRUE)
m[,interval_index_header] = cut(m[,opt$column], breaks=breaks, include.lowest=TRUE, right=TRUE, labels=FALSE)

# Format the interval
m[,interval_header] = sapply(m[,interval_header], formatInterval)

}


# Print output

output = ifelse(opt$output == "stdout", "", opt$output)
write.table(m, file=output, row.names=FALSE, quote=FALSE, col.names=opt$header, sep="\t")


# EXIT
quit(save='no')
