#!/usr/bin/env Rscript


options(stringsAsFactors=FALSE)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"),
	help="Matrix with expression values"),

make_option(c("-G", "--gene_pairs"), default="stdin",
	help="File or stdin. Two-column files with pairs of elements [default=%default]"),

make_option(c("-o", "--output"), default="cor.out.tsv",
	help="Output file name or stdout [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

make_option(c("-l", "--log"), action="store_true", default=FALSE,
	help="log10 of expression values [default=%default]"),

make_option(c("-p", "--pseudocount"), default=1e-04,
	help="pseudocount"),

make_option(c("-k", "--keep_na"), action="store_true", default=FALSE,
	help="do not replace NAs with 0 [default=%default]"),

make_option(c("-m", "--method"), default="pearson",
	help="method for correlation: pearson | spearman [default=%default]"),

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
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


# ======== #
# BEGIN    #
# ======== #

#input_files = strsplit(opt$input_files, ",")[[1]]
#labels = strsplit(opt$labels, ",")[[1]]
#vertical_lines = as.numeric(strsplit(opt$vertical_lines, ",")[[1]])


# Read data

if (opt$gene_pairs == "stdin") {
	pairs = read.table(file("stdin"), h=opt$header) 
} else {
	pairs = read.table(opt$gene_pairs, h=opt$header)
}

m = read.table(opt$input, h=T)

if (!(opt$keep_na)) {
	m = replace(m, is.na(m), 0)
}

if (opt$log) {
	m = log10(m + opt$pseudocount)
}


df = melt(as.matrix(m), value.name="value1", varnames=c("Var1", "variable"))
df = merge(df, pairs, by.x="Var1", by.y=colnames(pairs)[1])
df = merge(df, melt(as.matrix(m), value.name="value2", varnames=c("Var2", "variable")), 
	by.x=c(colnames(pairs)[2], "variable"), by.y=c("Var2","variable"))
df$split_by = apply(df[c(colnames(pairs)[2],"Var1")], 1, paste, collapse="_")
c = sapply(split(df, df$split_by), function(x) cor(x$value1, x$value2, u="p", m=opt$method))

df = merge(data.frame(c=c, split_by=names(c)), df)
df = unique(df[c("Var1", colnames(pairs)[2], "c")])

#print(head(df))
output = ifelse(opt$output == "stdout", "", opt$output)
write.table(df, output, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

#warnings()


q(save='no')

