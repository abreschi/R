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

make_option(c("-o", "--output"), default="quantile.out.tsv",
	help="Output file name. <stdout> for printing on stdout [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

make_option(c("-c", "--column"), default=1,
	help="The column with the values to divide in quantiles [default=%default]"),

make_option(c("-q", "--quantiles"), default=4, type="integer",
	help="Number of quantiles [default=%default]"),

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

quantile_header = sprintf("quantile_%s", colnames(m)[opt$column])
m[,quantile_header] = cut_number(m[,opt$column], opt$quantiles)

quant_index_header = sprintf("quant_index_%s", colnames(m)[opt$column])
m[,quant_index_header] = as.numeric(cut_number(m[,opt$column], opt$quantile))

# Print output

output = ifelse(opt$output == "stdout", "", opt$output)
write.table(m, file=output, row.names=FALSE, quote=FALSE, col.names=opt$header, sep="\t")


# EXIT
quit(save='no')
