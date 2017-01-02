#!/usr/bin/env Rscript 

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin", 
	help="Adjacency matrix in R format: the header is N-1 columns [default=%default]
		
		Example:
		
		sample1\tsample2\tsample3
		sample1\t1\t0\t1
		sample2\t1\t1\t0
		sample3\t0\t0\t1
	"),

make_option(c("--diag"), default=FALSE, action="store_true",
	help="Consider self-edges. [default=%default]"),

#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="Log-transform the data [default=%default]"),

#make_option(c("-p", "--pseudocount"), type="double", help="Pseudocount for the log [default=%default]"),

#make_option(c("--d1"), help="Design for voom"),
#make_option(c("-d", "--design"), help="Design for removing the batch effect (not including the batch effect)"),
make_option(c("-C", "--community"), help="Column in the metadata with the information of the community"),

make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-G", "--merge_mdata_on"), default="labExpId",
	help="Column in the metadata with the header of the input matrix [default=%default]"),

make_option(c("-o", "--output"), default="stdout", help="output file name [default=%default]"),
make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="verbose output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list, description="\nNormalize a matrix")
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}



# LIBRARIES

suppressPackageStartupMessages(library(igraph))
#suppressPackageStartupMessages(library(limma))


##--------------------##
## BEGIN              ##
##--------------------##


# read the matrix from the command line
if(opt$input_matrix == "stdin"){inF=file("stdin")}else{inF=opt$input_matrix}
m = read.table(inF, h=T, sep="\t", quote=NULL)

# Replace missing values with 0
m = replace(m, is.na(m), 0)

#if (opt$log) {m = log10(m + opt$pseudocount)}

g = graph.adjacency(as.matrix(m), mode='directed', diag=opt$diag, add.colnames=TRUE)
V = vertex.attributes(g)[[1]]

# =========================== Metadata =======================

merge_mdata_on = opt$merge_mdata_on
# read the metadata
mdata = read.table(opt$metadata, h=T, sep="\t", quote=NULL, comment.char="")
# Get the fields from the formula 
if (is.null(opt$community)) {
	cat("ERROR: please specify the batch variable\n")
	q(save='no')
}
fields = opt$community
mdata[opt$merge_mdata_on] <- gsub("[-,+]", ".", mdata[,opt$merge_mdata_on])
# Format the metadata
mdata = unique(mdata[unique(c(merge_mdata_on, fields))])
rownames(mdata) <- mdata[,merge_mdata_on]
mdata <- mdata[match(colnames(m), mdata[,merge_mdata_on]),, drop=FALSE]

community = mdata[match(colnames(m), mdata[, opt$merge_mdata_on]), opt$community]


modul = modularity(g, as.numeric(factor(community)))

# =================== OUTPUT ======================

out = paste(ecount(g), modul, sep="\t")
cat(out,"\n")
#cat(modul, "\n")

q(save='no')
