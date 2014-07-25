#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

##################
# OPTION PARSING #
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin [default=%default]"),

make_option(c("-o", "--output"), default="cutree.tsv",
	help="Output file name. Can be stdout [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="The input matrix has a header [default=%default]"),

make_option(c("-r", "--replace_na"), default=FALSE, action="store_true",
	help="Replace NAs with 0, before adding the pseudocount and applying the log if asked [default=%default]"),

make_option(c("-l", "--log10"), action="store_true", default=FALSE,
	help="Apply the log10 to the whole matrix as pre-processing step [default=%default]"),	

make_option(c("-p", "--pseudocount"), default=0.001,
	help="Pseudocount to add when applying the log [default=%default]"),

make_option(c("-k", "--nb_clusters"), 
	help="Number of desired clusters [default=%default]"),

make_option(c("-H", "--cluster_height"), 
	help="Height for cutting the tree [default=%default]"),

make_option(c("-d", "--dist"), default="euclidean",
	help="Distance measure for clustering. Supports spearman (s) and pearson (p) as distance metrics [default=%default]"),	

make_option(c("-c", "--hclust"), default="complete",
	help="Algorithm for the hierarchical clustering [default=%default]"),

make_option(c("--margin"), default=1,
	help="Cluster the rows (1) or the columns (2). Only rows can be clustered [default=%default]"),

#make_option(c("-B", "--iterations"), default=50,
#	help="Number of initializations to determine the best clustering [default=%default]"),
#
make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list = option_list,
	description = "Given a matrix cluster the element and cut the tree"
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


# TODO: cluster the columns

#suppressPackageStartupMessages(library("ggplot2"))


##############
# BEGIN
##############


if (opt$input == "stdin") {
	m = read.table(file("stdin"), h=T)
} else {
	m = read.table(opt$input, h=T)
}


if (opt$replace_na) {
	m = replace(m, is.na(m), 0)
}


if (opt$log10) {
	m = log10(m + opt$pseudocount)
}



# -------------------------- Distance -------------------------------

# Compute the distance between columns
if (opt$margin == 2) {
	if (opt$dist == "p" || opt$dist =="s") {
	    Dist = as.dist(1-cor(m, method=opt$dist, use="p"))
	} else {
	    Dist = dist(t(m), method=opt$dist)
	}
}

# Compute the distance between rows
if (opt$margin == 1) {
	if (opt$dist == "p" || opt$dist =="s") {
	    Dist = as.dist(1-cor(t(m), method=opt$dist, use="p"))
	} else {
	    Dist = dist(m, method=opt$dist)
	}
}


# -------------------------- Clustering ------------------------------

klust = hclust(Dist, method=opt$hclust)
K = cutree(klust, k=opt$nb_clusters, h=opt$cluster_height)
m$Kmeans = K



# OUTPUT

if (opt$output == "stdout") { 
	output = ""
} else {
	output = opt$output
}

write.table(m, output, quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t')

q(save='no')
