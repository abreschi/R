#!/usr/bin/env Rscript

##------------
## LIBRARIES
##------------

suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

option_list <- list(

make_option(c("-A", "--matrix_A"), 
	help="the matrix you want to subtract from, WITH header (A-B)"),

make_option(c("-B", "--matrix_B"), 
	help="the matrix you want to subtract, WITH header (A-B)"),

make_option(c("-r", "--replace_NA_with"), type="numeric",
	help="value you want to replace NA with, if null NAs are not replaced and difference will be NA"),

make_option(c("-o", "--output"), default="out.tsv",
	help="additional prefix for otuput [default=%default]"),

make_option(c("-e", "--expression"), 
	help="expression you want to corresponding cells of the matrices, e.g. A-B"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="verbose output [default=%default]")

)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if(opt$verbose) {print(opt)}

###################
# BEGIN           #
###################

A = read.table(opt$matrix_A, h=T)
B = read.table(opt$matrix_B, h=T)

if (!is.null(opt$replace_NA_with)) {
A <- replace(A, is.na(A), opt$replace_NA_with)
B <- replace(B, is.na(B), opt$replace_NA_with)
}

M = eval(parse(text=opt$expression))

write.table(M, opt$output, quote=F, row.names=T, sep="\t")

q(save='no')
