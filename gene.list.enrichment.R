#!/usr/bin/env Rscript


options(stringsAsFactors=FALSE)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin with a list of genes to test for enrichment [default=%default]"),

make_option(c("-d", "--db"), 
	help="A tab-separated database with the annotation. It has a header. First two columns are <gene_set> <gene_id>. \"%\" is the comment char."),

make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

#make_option(c("-U", "--all_genes_as_universe"), action="store_true", default=FALSE,
#	help="Use all genes in the list as universe, otherwise it uses only the genes annotated in db [default=%default]"),
#
make_option(c("-o", "--output"), default="gene.list.enrich.out.tsv",
	help="Output file name, stdout for printing on stdout [default=%default]"),

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
if (opt$verbose) {cat("DONE\n\n")}


# Wrapper for hypergeometric probaility

hyper.test = function(x) {
	x = as.numeric(x)
	q = x[1];
	m = x[2];
	n = x[4]-m;
	k = x[3];
	phyper(q, m, n, k, lower.tail=FALSE)
}



# ======== #
# BEGIN    #
# ======== #


# Read data

if (opt$input == "stdin") {inF = file("stdin")} else {inF = opt$input}
m = read.table(inF, h=opt$header)


db = read.table(opt$db, h=T, sep="\t", quote="", comment.char="%")

formula = as.formula(paste(colnames(db)[2],"~",colnames(db)[1]))
db_counts = setNames(aggregate(formula, db, length), c("feature", "total"))
list_db = db[db[,2] %in% m[,1], ]
list_counts = setNames(aggregate(formula, list_db, length), c("feature", "counts"))

df = merge(list_counts, db_counts, by="feature")
df$gene_list = length(unique(list_db[,2]))
df$universe = length(unique(db[,2]))
df$p.value = apply(df[-1], 1, hyper.test)
df$FDR = p.adjust(df$p.value, method="BH")

df = df[order(df$FDR),]

# Format and write output

df$p.value <- format(df$p.value, digits=2)
df$FDR <- format(df$FDR, digits=2)

output = ifelse(opt$output == "stdout", "", opt$output)
write.table(df, output, quote=FALSE, row.names=FALSE, sep="\t")


q(save='no')
