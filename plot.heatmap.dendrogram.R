#!/usr/bin/env Rscript


# This script is useful for:
# - clustering the samples by gene expression given a matrix
###########
# VARIABLES
###########

options(stringsAsFactors=F)
hcluster <- "complete"
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
cbbPalette <- rgb(matrix(c(0,0,0,0,73,73,0,146,146,255,109,182,255,182,119,73,0,146,0,109,219,182,109,255,109,182,255,182,219,255,146,0,0,146,73,0,219,209,0,36,255,36,255,255,109), ncol=3, byrow=T), max=255)

my_palette <- colors()[c(12,17,24,26,31,33,35,50,51,56,62,84,92,116,139,153,371,490,645)]

pseudocount <- 1e-04


##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-m", "--metadata"), help="one or more tsv file(s) with metadata on matrix experiment"),
make_option(c("-d", "--dist"), help="distance measure between samples. Choose among <p> (pearson), <s> (spearman). [default=%default]", default="p"),
make_option(c("-c", "--hclust"), help=sprintf("hierarchical clustering method. Choose among <complete>, <average> [default=%s]",hcluster), default=hcluster),
make_option(c("-C", "--color_by"), help="metadata by which colouring the variables [default=NA]", type='character', default=NA ),
make_option(c("-o", "--output_suffix"), help="additional suffix for the output [default=%default]", default="out"),
make_option(c("-f", "--fields"), help="write the header of the fields you want in the labels, comma-separated [default=labExpId]", default="labExpId")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

print(opt)

##------------
## LIBRARIES
##------------ 

cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
library(RColorBrewer)
cat("DONE\n\n")


##--------------------##
## CLUSTERING SAMPLES ##
##--------------------##

# read the matrix from the command line
m = read.table(opt$input_matrix, h=T)

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}


# substitute the matrix with its log if required by the user
if (opt$log) {
m = log2(replace(m, m==0 | is.na(m), opt$pseudocount))}


# read metadata, one or more table to be merged on cell
metadata = strsplit(opt$metadata, ",")[[1]]
for (i in seq_along(metadata)) {
    mdata = read.table(metadata[i], h=T, sep="\t", row.names=NULL);
	if ('labExpId' %in% colnames(mdata)) {
    mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))}
#	mdata$labExpId <- sapply(mdata$labExpId, function(x) if(any(startsWith(x,0:9))){return(sprintf('X%s', x))})
	if ( i==1 ) {
        new_mdata = mdata
    }else{
    new_mdata = merge(mdata, new_mdata, by=c('cell'))
    }
}




# compute the correlation according to user's choice
if (opt$dist == 'p') {
m_cor = cor(m, use='p', method='p')
m_dist = 1 - abs(m_cor)
}

if (opt$dist == 's') {
m_cor = cor(m, use='p', method='s')
m_dist = 1 - abs(m_cor)
}


# concatenate the labels from the metadata
labels = sapply(colnames(m_cor), function(x) paste(unique(subset(new_mdata, labExpId==x)[strsplit(opt$fields, ",")[[1]]]), collapse='_'))


# dendrogram
hcluster = opt$hclust
colLab = function(n,sample_ids=NULL, colors=NULL,new_sample_names=NULL) {
  if(is.leaf(n)) {
    a <- attributes(n);
	new_lab = new_sample_names[which(sample_ids == a$label)]
    node_col = colors[which(sample_ids == a$label)]
    attr(n,'nodePar') <- list(lab.col=node_col, lab.cex= 1/log10(length(sample_ids)))
	attr(n,'label') <- new_lab
    };n}


# passing the colors
color_by = strsplit(opt$color_by, ',')[[1]]
if (is.na(opt$color_by)) {rc = rep('black', ncol(m_cor))}else{
rc = cbbPalette[as.numeric(as.factor(sapply(colnames(m_cor), 
function(x) paste(unique(subset(new_mdata, labExpId==x)[,color_by]), collapse=','))))];
legend = sapply(colnames(m_cor), function(x) paste(unique(subset(new_mdata, labExpId==x)[,color_by]), collapse=','))}





# print output file name
output = sprintf('heatmap.dendrogram.log_%s.psd_%s.dist_%s.hclust_%s.%s.pdf',
opt$log, opt$pseudocount, opt$dist, opt$hclust, opt$output_suffix)

#pdf(output, width=ncol(m_cor)/10+5, height=nrow(m_cor)/10+10)
right_omi = log2(max(nchar(labels)))/2+0.5
pdf(output, width=7+right_omi, height=7)
par(omi=c(0.2,0.2,0.2,right_omi), xpd=NA, mai=c(0,0,0,0), fin=c(1.5,1.5), cex.main=0.75)
heatmap.2(m_cor, distfun=function(x) dist(x,method='euclidean'),
hclustfun=function(x) hclust(x,method=hcluster),trace='none',
cexRow = 1/log10(nrow(m_cor)), cexCol = 1/log10(ncol(m_cor)),
col=rev(heat.colors(16)),keysize=1,labCol='',
#main=opt$dist,
main="",
labRow=labels, RowSideColors=rc
#lhei=c(1,3)
)
legend(x=1, y=1, legend=unique(legend), fill=unique(rc))

#par(mar=c(7,2,2, right_oma),cex.sub=2,oma=c(5,1,1,1))
par(mai=c(1,0.2,0.2,right_omi*2), cex.sub=2, omi=c(0.4,0.1,0.1,0.1))
dend = as.dendrogram(hclust(as.dist(m_dist),hcluster))
plot(dendrapply(dend, colLab, sample_ids=colnames(m_cor), colors=rc, new_sample_names=labels), horiz=T)
#title(sub=sprintf('dist=%s,hclust=%s',opt$dist,hcluster))

dev.off()

q(save='no')


