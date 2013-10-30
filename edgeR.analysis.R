#!/usr/bin/env Rscript 

# -- Variables --

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


#########################
# Default OPT parameters
########################

opt = list()
opt$n_samples = 2
opt$cpm = 1

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix with READ COUNTS you want to analyze"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-f", "--fields"), help="choose the fields you want to use in the differential expression, comma-separated"),
make_option(c("-F", "--formula"), help="a formula for the differential expression [default=%default]", default="~."),
make_option(c("-C", "--cpm"), help="threshold for cpm [default=%default]", type="integer", default=1),
make_option(c("-s", "--n_samples"), help="minimum number of samples with cpm > \"cpm\" [default=%default]", type="integer", default=2),
make_option(c("-c", "--coefficient"), help="the coefficient for the comparison [default=last field]", type='integer'),
make_option(c("-o", "--output"), help="additional suffix for output [default=%default]", default='out'),
make_option(c("-d", "--output_dir"), help="choose the output directory [default=%default]", default="./"),
make_option(c("-g", "--genes"), help='a file with a list of genes to filter', type='character', default=NA)
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)


##------------
## LIBRARIES
##------------ 


cat('Libraries loading... ')

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library('edgeR'))

cat('DONE\n\n')



##-------------------------##
## Differential Expression ##
##-------------------------##

plotMA = function(x, fdr) {plot(sort(x$logCPM), x$logFC[order(x$logCPM)]);
points(subset(x, FDR<fdr)[c("logCPM","logFC")], col='red')}


# read the matrix from the command line
m = read.table(opt$input_matrix, h=T)
genes = rownames(m)
m = (apply(m, 2, as.integer))
rownames(m) <- genes

# read the matrix
if (!is.na(opt$genes)) {
fgenes = as.character(read.table(opt$genes, h=F)$V1)
m = m[intersect(rownames(m),fgenes),]
}

# read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')

# specify the design to the program
fields = strsplit(opt$fields, ",")[[1]]

# ONLY ONE CONDITION

if (length(fields) == 1) {
condition = as.factor(sapply(colnames(m), function(x) unique(subset(mdata, labExpId == x)[,fields])))

# create count object for edgeR
M = DGEList(counts=na.omit(m), group = condition)

# normalisation
cpm.m <- cpm(M)
M <- M[rowSums(cpm.m > opt$cpm) >= opt$n_samples,]
M$samples$lib.size <- colSums(M$counts)
M <- calcNormFactors(M)

# variance estimation
M <- estimateCommonDisp(M, verbose=T)
M <- estimateTagwiseDisp(M)

# calling differential expression
et <- exactTest(M)
res = topTags(et, n=nrow(et))$table

# MULTIPLE CONDITIONS

}else{
fields = strsplit(opt$fields, ",")[[1]]
design_df = unique(subset(mdata, labExpId %in% colnames(m))[c("labExpId",fields)])
design_df = design_df[order(design_df$labExpId),]
design = model.matrix(as.formula(opt$formula), design_df[fields])
print(colnames(design))
rownames(design) <- design_df$labExpId

cat('\n')


# create count object for edgeR
M = DGEList(counts=na.omit(m))

# normalisation
cpm.m <- cpm(M)
M <- M[rowSums(cpm.m > opt$cpm) >= opt$n_samples,]
M$samples$lib.size <- colSums(M$counts)
M <- calcNormFactors(M)

# variance estimation
M <- estimateGLMCommonDisp(M, design, verbose=T)
M <- estimateGLMTrendedDisp(M, design)
M <- estimateGLMTagwiseDisp(M, design)


# calling differential expression
fit <- glmFit(M, design)
coeff <- ifelse(is.null(opt$coefficient), ncol(fit), opt$coefficient)
lrt <- glmLRT(fit, coef=coeff)
res = topTags(lrt, n=nrow(lrt))$table
}

# otuput results to a file
output = sprintf("%s/edgeR.cpm%s.n%s.%s", opt$output_dir, opt$cpm, opt$n_samples, opt$output)
write.table(res, file=sprintf("%s.tsv", output), quote=F, sep='\t',row.names=T)
#write.table(cpm.m, file=sprintf("cpm.%s.tsv",opt$output), quote=F, sep='\t',row.names=T)

#------#
# PLOT #
#------#

pdf(sprintf('%s.pdf', output))
# MA plot
plotMA(res, 0.01)
# plot dispersion estimates
plotBCV(M)
dev.off()

q(save='no')
