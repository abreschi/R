#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

##################
# OPTION PARSING #
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin [default=%default]"),

make_option(c("-o", "--output"), default="anova.out.tsv",
	help="Output file name. Can be stdout [default=%default]"),

make_option(c("-m", "--metadata"),
	help="Matrix with the metadata"),

make_option(c("--merge_mdata_on"), default="labExpId",
	help="Metadata field which contains the column names of the input matrix [default=%default]"),

make_option(c("-F", "--factors"), 
	help="Factors for anova, can also be interactions, e.g. value~cell+organism"),

make_option(c("--p_adj"), default="BH",
	help="Method for correcting the pvalue for multiple testing [default=%default]"),

make_option(c("-l", "--log10"), action="store_true", default=FALSE,
	help="Apply the log10 to the whole matrix as pre-processing step [default=%default]"),	

make_option(c("-p", "--pseudocount"), default=0.001,
	help="Pseudocount to add when applying the log [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))


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


# Read the metadata
mdata = read.table(opt$metadata, h=T, sep="\t")
mdata[,opt$merge_mdata_on] = gsub(",", ".", mdata[,opt$merge_mdata_on])
mdata_col = unique(c(opt$merge_mdata_on, strsplit(opt$factors, "[+*:]")[[1]]))
mdata = unique(mdata[,mdata_col])

# Read the formula
F = as.formula(sprintf("value~%s", opt$factors))

#m = m[1:100,]

# Apply anova on each gene
res = t(sapply(1:nrow(m), 
	function(i) {
		mm = suppressWarnings(melt(m[i,]))
		tmp = merge(mm, mdata, by.x="variable", by.y=opt$merge_mdata_on);
		aov_res = anova(lm(F, tmp));
		return(unlist(aov_res[,-1]))
	}
))

tmp = merge(melt(m[1,]), mdata, by.x="variable", by.y=opt$merge_mdata_on);
aov_res = anova(lm(F, tmp));
labels = apply(expand.grid(rownames(aov_res), c("SS", "MeanSq", "F", "pvalue")), 1, paste, collapse="_") 

res = data.frame(res)
colnames(res) <- labels
rownames(res) <- rownames(m)


# Adjust pvalue
for (i in grep("pvalue", colnames(res))) {
	adj_header = paste(colnames(res)[i], "adj", sep=".")
	res[,adj_header] = p.adjust(res[,i], method=opt$p_adj)
}
	
	

# OUTPUT

if (opt$output == "stdout") { 
	output = ""
} else {
	output = opt$output
}

write.table(res, output, quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t')

q(save='no')
