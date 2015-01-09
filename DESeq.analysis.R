#!/usr/bin/env Rscript 


##------------
## LIBRARIES
##------------ 

cat('Loading libraries... ')

#suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
#suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library('DESeq'))

cat('DONE\n\n')

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix with READ COUNTS you want to analyze"),
make_option(c("-r", "--replace_NAs"), action="store_true", default=FALSE, help="replace NAs with 0 [default=%default]"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-f", "--fields"), help="choose the fields you want to use in the differential expression, comma-separated"),
make_option(c("-o", "--output"), help="output file name"),
make_option(c("-v", "--verbose"), action="store_true", default=FALSE)
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if(opt$verbose) {print(opt)}



##--------------------##
## CLUSTERING SAMPLES ##
##--------------------##

# read the matrix from the command line
m = read.table(opt$input_matrix, h=T)
genes = rownames(m)
m = (apply(m, 2, as.integer))
rownames(m) <- genes

if (opt$replace_NAs) {m = replace(m, is.na(m), 0)}

# read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')
mdata["labExpId"] <- gsub("[,-]", ".", mdata[,"labExpId"])


# specify the design to the program
fields = strsplit(opt$fields, ",")[[1]]
#mdata = unique(mdata[c("labExpId", opt$fields)])
#print(match(colnames(m), mdata[,"labExpId"]))
if (length(fields) == 1) {
condition = factor(sapply(colnames(m), function(x) unique(subset(mdata, labExpId == x)[,opt$fields])))
}else{print('cannot handle multiple fields yet');q(save='no')}


#colData = mdata[mdata[,"labExpId"] %in% colnames(m) & mdata[,"view"]=="Alignments",]
#rownames(colData) <- colData[,"labExpId"]
#colData = colData[match(colnames(m), rownames(colData)),]
#dds = DESeqDataSetFromMatrix(countData = m, colData= colData, design=~cell)



# create count object for DESeq
cds = newCountDataSet(m, condition)

# normalization
if(opt$verbose) {cat("Estimating size factors... ")}
cds = estimateSizeFactors(cds)
if(opt$verbose) {cat("DONE\n")}

# variance estimation
if(opt$verbose) {cat("Estimating dispersions... ")}
cds = estimateDispersions(cds)
if(opt$verbose) {cat("DONE\n")}

# plot dispersion estimates
plotDispEsts = function(x) {
norm_means = rowMeans(counts(x, normalized=T))
plot(sort(norm_means), fitInfo(x)$perGeneDispEsts[order(norm_means)], log='xy')
lines(sort(norm_means), fData(x)$disp_pooled[order(norm_means)],  col=2)
}

# calling differential expression
if(opt$verbose) {cat("Binomial test... ")}
res = nbinomTest(cds, levels(condition)[1], levels(condition)[2])
if(opt$verbose) {cat("DONE\n")}

# plot MA
plotMA = function(x) {
plot(sort(log10(x$baseMean)), x$log2FoldChange[order(x$baseMean)])
points(log10(subset(x, padj<=0.01)[,'baseMean']), subset(x, padj<=0.01)[,'log2FoldChange'], col='red')
}



# output on file
toTable = data.frame("id"=res[,1], "log2FC"=res[,6], "log2CPM"=log2(res[,2]), "PValue"=res[,7], "FDR"=res[,8])
write.table(toTable, file=opt$output, quote=F, sep='\t', row.names=F)

q(save='no')


# for debugging
opt = list()
opt$input_matrix = '~/Documents/blueprint/pilot/Flux/bp.human.gene.reads.idr_01.thr_0.names_False.tsv'
opt$metadata = '~/Documents/blueprint/pilot/bp_rna_dashboard_temp.columns'
opt$fields = 'cell'
opt$genes = '~/Documents/blueprint/pilot/gen15.gene.pc.txt'

