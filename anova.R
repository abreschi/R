#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

##################
# OPTION PARSING #
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="Input matrix (R-friendly, the header has n-1 columns [default=%default]"),

make_option(c("-o", "--output"), default="stdout",
	help="Output file name. [default=%default]"),

make_option(c("-l", "--log10"), action="store_true", default=FALSE,
	help="Apply the log10 to the whole matrix as pre-processing step [default=%default]"),	

make_option(c("-p", "--pseudocount"), default=0.01,
	help="Pseudocount to add when applying the log [default=%default]"),

make_option(c("-r", "--replace_NA"), default=FALSE, action="store_true",
	help="Replace NA with 0 [default=%default]"),

make_option(c("-m", "--metadata"),
	help="Matrix with the metadata. It contains the information about the columns of the input matrix"),

make_option(c("--merge_mdata_on"), default="labExpId",
	help="Metadata field which contains the column names of the input matrix [default=%default]"),

make_option(c("-F", "--factors"), 
	help="Factors for anova (right part of the formula after \"~\"), can also be interactions, e.g. value~cell+organism"),

make_option(c("--p_adj"), default="BH",
	help="Method for correcting the pvalue for multiple testing [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list,
	description = "
	Computes anova for each row of the input matrix given a design model. 
	The SS for each factor of the design can be used to compute the 
	proportion of variance explained by that factor

	Requires package \"reshape2\"
	"
)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


suppressPackageStartupMessages(library("reshape2"))


##############
# BEGIN
##############

# Read input
#input = ifelse(opt$input == "stdin", file("stdin"), opt$input) # why is this line never working?
if (opt$input == "stdin") {
        input = file("stdin")
} else {
        input = opt$input
}
m = read.table(input, h=T)
if (opt$verbose) {
	cat("Data sample:\n")
	print(head(m))
}

# Process the data according to the user
if (opt$replace_NA) {m = replace(m, is.na(m), 0)}
if (opt$log10) {m = log10(m + opt$pseudocount)}


# Read the metadata
mdata = read.table(opt$metadata, h=T, sep="\t", quote=NULL, check.names=F)
mdata[,opt$merge_mdata_on] = gsub("[,:-]", ".", mdata[,opt$merge_mdata_on])
mdata_col = unique(c(opt$merge_mdata_on, strsplit(opt$factors, "[+*:]")[[1]]))
mdata = unique(mdata[,mdata_col])
if (opt$verbose) {
	cat("Metadata sample:\n")
	print(head(mdata))
}

# Read the formula
F = as.formula(sprintf("value~`%s`", opt$factors))

#m = m[1:100,]

# Apply anova on each gene
res = t(sapply(1:nrow(m), 
	function(i) {
		mm = suppressMessages(melt(m[i,]))
		tmp = merge(mm, mdata, by.x="variable", by.y=opt$merge_mdata_on);
#		print(tmp)
		aov_res = anova(lm(F, tmp));
		return(unlist(aov_res[,-1]))
	}
))

tmp = merge(melt(m[1,]), mdata, by.x="variable", by.y=opt$merge_mdata_on);
aov_res = anova(lm(F, tmp));
labels = apply(expand.grid(rownames(aov_res), c("SS", "MeanSq", "F", "pvalue")), 1, paste, collapse="_") 

res = data.frame(res)
colnames(res) <- labels


# Adjust pvalue
for (i in grep("pvalue", colnames(res))) {
	adj_header = paste(colnames(res)[i], "adj", sep=".")
	res[,adj_header] = p.adjust(res[,i], method=opt$p_adj)
}

res = sapply(res, round, 4)
rownames(res) <- rownames(m)

# If the variance is zero set the results to NA
res[apply(m, 1, var)==0,] <- NA


# OUTPUT

output = ifelse(opt$output == "stdout", "", opt$output)
write.table(res, output, quote=FALSE, col.names=TRUE, row.names=TRUE, sep='\t')

q(save='no')
