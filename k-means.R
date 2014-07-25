#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

##################
# OPTION PARSING #
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin [default=%default]"),

make_option(c("-o", "--output"), default="Kmeans.tsv",
	help="Output file name. Can be stdout [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="The input matrix has a header [default=%default]"),

make_option(c("-l", "--log10"), action="store_true", default=FALSE,
	help="Apply the log10 to the whole matrix as pre-processing step [default=%default]"),	

make_option(c("-p", "--pseudocount"), default=0.001,
	help="Pseudocount to add when applying the log [default=%default]"),

make_option(c("-k", "--nb_clusters"), default=3,
	help="Number of desired clusters [default=%default]"),

make_option(c("-B", "--iterations"), default=50,
	help="Number of initializations to determine the best clustering [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


suppressPackageStartupMessages(library("ggplot2"))


##############
# BEGIN
##############


if (opt$input == "stdin") {
	m = read.table(file("stdin"), h=T)
} else {
	m = read.table(opt$input, h=T)
}

if (opt$log10) {
	m = log10(m + opt$pseudocount)
}


set.seed(123)

# Find the clusters from multiple random initializations
Klist = replicate(opt$iterations, kmeans(m, opt$nb_clusters), simplify=F)

# Choose the best
K = Klist[[which.max(sapply(1:length(Klist), function(i) {Klist[[i]]$betweenss/Klist[[i]]$totss}))]]$cluster

m$Kmeans = K


# OUTPUT

if (opt$output == "stdout") { 
	output = ""
} else {
	output = opt$output
}

write.table(m, output, quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t')

q(save='no')
