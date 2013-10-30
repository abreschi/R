#!/usr/bin/env Rscript

# == NOTES ==
# 1) the gtf should NOT have the trailing ";"
# 2) the list of files with the exon read counts should have full paths
# 3) the order of the id should be the same as the files


##------------
## LIBRARIES
##------------ 

cat('Loading libraries... ')

#suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library('DEXSeq'))

cat('DONE\n\n')

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

# =====================
# == DEBUG OPTIONS ==
# ====================

opt = list()
opt$metadata = "~/Documents/blueprint/pilot/bp_rna_dashboard_mmaps.crg.tsv"
opt$countfiles = "ERR180942.exon.count.txt,ERR180943.exon.count.txt,ERR180944.exon.count.txt,ERR180945.exon.count.txt,ERR180948.exon.count.txt,ERR180950.exon.count.txt,ERR180951.exon.count.txt,ERR186015.exon.count.txt,ERR230581.exon.count.txt,ERR232403.exon.count.txt,ERR232404.exon.count.txt,ERR244135.exon.count.txt"
opt$labExpId = "ERR180942,ERR180943,ERR180944,ERR180945,ERR180948,ERR180950,ERR180951,ERR186015,ERR230581,ERR232403,ERR232404,ERR244135"
opt$DE_factor = "cell"
opt$annotationfile = "~/Documents/db/human/gencode15/Long/DEXSeq/gen15.long.DEXSeq.gtf"


##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--countfiles"), help="list of files with the counts for exons"),
#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
#make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-I", "--labExpId"), help="list of ids in the same order as the files"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-f", "--DE_factor"), help="choose the field you want to use for the DE"),
make_option(c("-a", "--annotationfile"), help="gtf file"),
make_option(c("-o", "--output"), help="output file name"),
#make_option(c("-f", "--fill_by"), help="choose the color you want to fill by [default=NA]", type='character', default=NA)
make_option(c("-g", "--genes"), help='a file with a list of genes to filter', type='character', default=NA)
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)


# == BEGIN ==

minCount = 24
maxExon = 20
nCores = 6
progressFile = "DEXSeq.progress"

mdata = read.table(opt$metadata, h=T, stringsAsFactors=T)
countfiles = strsplit(opt$countfiles, ",")[[1]]
labExpId = strsplit(opt$labExpId, ",")[[1]]

# Prepare the design matrix from metadata
df = unique(subset(mdata, labExpId %in% labExpId, select=c("labExpId", opt$DE_factor)))
df = df[match(labExpId, df$labExpId),]
rownames(df) <- df$labExpId
colnames(df)[colnames(df) == opt$DE_factor] <- "condition"

ecs = read.HTSeqCounts(countfiles = countfiles, design = df, flattenedfile=opt$annotationfile)
sampleNames(ecs) <- labExpId

# Run the test on a subset of genes
test = subsetByGenes(ecs, genes=sample(geneIDs(ecs), 500))

# Estimate the scaling factors
test <- estimateSizeFactors(test)

# Estimate dispersions
test <- estimateDispersions(test, minCount=minCount, maxExon=maxExon, nCores=nCores, file=progressFile)

# Fit dispersions
test <- fitDispersionFunction(test)

# Make sure the regression line is ok
plot(log10(rowMeans(counts(test, normalized=T))), log10(fData(test)$dispBeforeSharing))
points(log10(rowMeans(counts(test, normalized=T))), log10(fData(test)$dispFitted), col='red')

# Test the DE
test <- testForDEU(test, nCores=nCores)

# Store the results in a data.frame
res = DEUresultTable(test)
plotDEXSeq(test, "ENSG00000139505.10", displayTranscripts=T, cex.axis=1.2, cex=1.3, lwd=2, legend=T, norCounts=T)


