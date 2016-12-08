#!/usr/bin/env Rscript

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), type="character", default='stdin',
	help="tab-separated file. Can be stdin [default=%default]"),

make_option(c("-o", "--output"), default="stdout", 
	help="output file name with extension [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE, 
	help="The file has header [default=%default]"),

make_option(c("-V", "--values"), type='numeric', default=1,
	help="Column index with the values [default=%default]"),

make_option(c("-f", "--factor"), type='numeric', default=1,
	help="Column index with the factor [default=%default]"),

make_option(c("-v", "--verbose"), action='store_true', default=FALSE,
	help="Verbose output [default=%default]")

#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list = option_list,
	description = "Compute wilcox.test for each pair of group in the column factor"
)


arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


##------------
## LIBRARIES
##------------
 
#if (opt$verbose) {cat("Loading libraries... ")}
#suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(ggplot2))
##suppressPackageStartupMessages(library(plyr))
#if (opt$verbose) {cat("DONE\n\n")}


# ========
# BEGIN
# ========


# Read data
if (opt$input == "stdin") {input=file('stdin')} else {input=opt$input}
m = read.table(input, h=opt$header, sep="\t")
if(opt$verbose) {print(head(m))}

# Read the axes
values = colnames(m)[opt$values]
fact = colnames(m)[opt$factor]
# Build the formula
form = as.formula(sprintf("%s~%s", values, fact))

if(opt$verbose) {
	cat("Formula: ")
	print(form)
}

# Vector of groups
groups = unique(as.character(m[,fact]))

# Perform mann-whitney on all pairs of groups
res = apply(combn(groups, 2), 2, 
	function(x) {
		df = m[m[,fact] %in% x,];
		res.test = wilcox.test(form, data=df);
		data.frame(
			group1 = x[1],
			group2 = x[2],
			W = res.test$statistic,
			p.value = format(res.test$p.value, digits=4)
		)
	}
)

# Concatenate the data.frame
res = do.call(rbind, res)

# OUTPUT
output = ifelse(opt$output == "stdout", "", opt$output)
write.table(res, file=output, row.names=FALSE, quote=FALSE, sep="\t")

# EXIT
quit(save='no')
