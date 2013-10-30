
# This script is useful for:
# - clustering the samples by gene expression given a matrix

##------------
## LIBRARIES
##------------ 
cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(RColorBrewer))
cat("DONE\n")

###########
# VARIABLES
###########

options(stringsAsFactors=F)
pseudocount <- 1e-04
hcluster <- "complete"
dist <- "euclidean"
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
custom_palette <- colors()[c(12,17,24,26,31,33,35,50,51,56,62,84,92,116,139,153,371,490,645)]

my_palette = cbbPalette


##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-s", "--scale_heatmap"), action="store_true", default=FALSE, help="apply scaling by row [default=FALSE]"),
make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-d", "--dist"), help=sprintf("distance measure between samples. Choose among <euclidean> [default=%s]",dist), default=dist),
make_option(c("-c", "--hclust"), help=sprintf("hierarchical clustering method. Choose among <complete>, <average> [default=%s]",hcluster), default=hcluster),
make_option(c("-G", "--genes"), help="one column of gene ids you want to plot the expression of, no header. [optional]"),
make_option(c("-cm", "--col_mdata"), help="tsv file with metadata related to the column header"),
make_option(c("-rm", "--row_mdata"), help="tsv file with metadata related to the row names. Column with the gene id must have header gene_id"),
make_option(c("-cc", "--col_color_by"), help="metadata by which colouring the columns [default=NA]", type='character', default=NA ),
make_option(c("-rc", "--row_color_by"), help="metadata by which colouring the rows [default=%default]", type='character', default=NA),
make_option(c("-cl", "--col_labels"), help="write the header of the fields you want in the column labels, comma-separated, NA for no labels [default=colnames]"),
make_option(c("-rl", "--row_labels"), help="the names of the metadata field you want in the row labels, comma-separated, NA for no labels"),
make_option(c("-D", "--dendro"), help="where to apply cluster. <both>, <none>, <row>, <col> [default=%default]", default="both"),
make_option(c("-o", "--output_suffix"), help="additional suffices for output before extension [default=out]", default='out'),
make_option(c("-t", "--transpose"), help="plot the values in transposed order", action="store_true", default=FALSE)
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

print(opt)


##--------------------##
## CLUSTERING SAMPLES ##
##--------------------##

# read the matrix from the command line
m = read.table(opt$input_matrix, h=T)

# substitute the matrix with its log if required by the user
if (opt$log) {
m = log10(replace(m, m==0 | is.na(m), opt$pseudocount))}




# ------------- SELECT GENES ------------------------------------------

# select the genes in the list provided by the user
if (!is.null(opt$genes)) { 
fgenes = read.table(opt$genes, col.names='gene', h=FALSE)
m <- m[intersect(unique(fgenes$gene), rownames(m)),]
}



# ------------- COLUMN METADATA --------------------------------------

# read the metadata for the columns from the metadata file if provided
if (!is.null(opt$col_mdata)) {
col_mdata = read.table(opt$col_mdata, h=T, sep='\t')
# substitute commas and minuses in the labExpId field
col_mdata$labExpId <- sapply(col_mdata$labExpId, function(x) gsub("[,-]", ".", x))
}

# concatenate the labels from the metadata if specified
if (!is.null(opt$col_labels)) {
if (opt$col_labels!="NA") {
col_labels = sapply(colnames(m), function(x) paste(unique(subset(col_mdata, labExpId==x)[strsplit(opt$col_labels, ",")[[1]]]), collapse='_'))
bottom_omi <- max(nchar(col_labels))*0.13   # strwidth("O",u="i")
}else{col_labels = ''; bottom_omi <- 0.5}
}else{col_labels = colnames(m); bottom_omi <- max(nchar(col_labels))*0.13}

# passing the column colors if specified
if (is.na(opt$col_color_by)) {cc = rep('black', ncol(m))}else{
cc = my_palette[as.numeric(as.factor(sapply(colnames(m), function(x) unique(subset(col_mdata, labExpId==x)[,opt$col_color_by]))))]}



# ------------ ROW METADATA ----------------------------------------

# read the metadata for the rows from the metadata file if provided
if (!is.null(opt$row_mdata)) {
row_mdata = read.table(opt$row_mdata, h=T, sep='\t')
# substitute commas and minuses in the labExpId field
row_mdata$gene_id <- sapply(row_mdata$gene_id, function(x) gsub("[,-]", ".", x))
}

# concatenate the labels from the metadata if specified
if (!is.null(opt$row_labels)) {
if (opt$row_labels!="NA") {
row_labels = sapply(rownames(m), function(x) paste(unique(subset(row_mdata, gene_id==x)[strsplit(opt$row_labels, ",")[[1]]]), collapse='_'))
right_omi = max(nchar(row_labels))*0.13  # strwidth("O",u="i")
}else{row_labels = ''; right_omi <- 0.5}
}else{row_labels = rownames(m); right_omi = max(nchar(row_labels))*0.13}

# passing the column colors if specified
if (is.na(opt$row_color_by)) {rc = rep('black', nrow(m))}else{
N_COLORS = length(unique(row_mdata[,opt$row_color_by]))
rc = brewer.pal(9, "Paired")[as.numeric(as.factor(sapply(rownames(m), function(x) unique(subset(row_mdata, gene_id==x)[,opt$row_color_by]))))];
legend = sapply(rownames(m), function(x) paste(unique(subset(row_mdata, gene_id==x)[,opt$row_color_by]), collapse=','))}


# ----------------- DENDROGRAMS -----------
if (opt$dendro == "both") {Rowv = T; Colv = T}
if (opt$dendro == "none") {Rowv = F; Colv = F}
if (opt$dendro == "col") {Rowv = F; Colv = T; if (opt$transpose){Rowv=!Rowv; Colv=!Colv}}
if (opt$dendro == "row") {Rowv = T; Colv = F; if (opt$transpose){Rowv=!Rowv; Colv=!Colv}}


# -------------- OUTPUT ------------------------------
#
# print output file name
output = sprintf('heatmap.log_%s.pscnt_%s.dist_%s.hclust_%s.%s.pdf', 
opt$log, opt$pseudocount, opt$dist, opt$hclust, opt$output_suffix)



if (!(opt$transpose)) {
# interpret user choice about scaling
scale <- ifelse(opt$scale_heatmap, 'row', 'none')
#pdf(output, width=15, height=20)
# genes are rows and samples are columns
# strheight("M\n", u="i") --> 0.32
pdf(output, width=ncol(m)*0.22+right_omi+3, height=nrow(m)*0.22+bottom_omi)
par(omi=c(bottom_omi,0.2,0.2,right_omi), xpd=NA, mai=c(0,0,0,0))
heatmap.2(as.matrix(m), 
Rowv = Rowv,
Colv = Colv,
scale=scale,
na.color='grey',
distfun=function(x) dist(x,method=opt$dist),
hclustfun=function(x) hclust(x,method=hcluster),
trace='none',
col=redgreen(75),
keysize=1,
#main=opt$dist,
main="",
labCol=col_labels, 
labRow=row_labels,
ColSideColors=cc,
RowSideColors=rc,
cexCol=1,
cexRow=1
)
if (!is.na(opt$row_color_by)) {
legend(x=1, y=1, legend=unique(legend), fill=unique(rc), cex=0.75)}
#box("plot")
#box("figure", col="green")
#box("outer", col="blue")
dev.off()
}

if (opt$transpose) {
# interpret user choice about scaling
scale <- ifelse(opt$scale_heatmap, 'column', 'none')
pdf(output, width=15+bottom_omi, height=7+right_omi)
par(omi=c(right_omi,0.2,0.2,bottom_omi), xpd=NA, mai=c(0,0,0,0), fin=c(1.5,1.5))
#pdf(output, width=15, height=10)
# genes are columns and samples are rows
#par(oma=c(1.5, 2, 2, bottom_oma), xpd=NA)
heatmap.2(t(as.matrix(m)), 
Rowv = Rowv,
Colv = Colv,
scale=scale,
na.color='grey',
distfun=function(x) dist(x,method=opt$dist),
hclustfun=function(x) hclust(x,method=hcluster),trace='none',
cexRow = 0.2 + 1/log10(ncol(m)),
cexCol = 0.2 + 1/log10(nrow(m)),
col=redgreen(75),
keysize=1,
main="",
labRow=col_labels, 
labCol=row_labels,
RowSideColors=cc, 
ColSideColors=rc,
lhei=c(ncol(m)*1/5,ncol(m)*4/5))
if (!is.na(opt$row_color_by)) {
legend(x=0, y=0, legend=legend, fill=unique(rc))}
dev.off()
}


q(save='no')






## dendrogram
#hcluster = opt$hclust
#colLab = function(n,sample_ids=NULL, colors=NULL,new_sample_names=NULL) {
#  if(is.leaf(n)) {
#    a <- attributes(n);
#	new_lab = new_sample_names[which(sample_ids == a$label)]
#    node_col = colors[which(sample_ids == a$label)]
#    attr(n,'nodePar') <- list(lab.col=node_col, lab.cex=1+4/5)
#	attr(n,'label') <- new_lab
#    };n}


par(mar=c(7,2,2,50),cex.sub=2,oma=c(5,1,1,1))
dend = as.dendrogram(hclust(as.dist(m_dist),hcluster))
plot(dendrapply(dend, colLab, sample_ids=colnames(m_cor), colors=rc, new_sample_names=labels), horiz=T)
title(sub=sprintf('dist=%s,hclust=%s',opt$dist,hcluster))

dev.off()

q(save='no')


