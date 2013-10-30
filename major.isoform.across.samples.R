
#############################
# LIBRARIES
#####################

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix with transcript expression you want to analyze"),
make_option(c("-a", "--annotation"), help="two-column file with gene and tx ids"),
make_option(c("-o", "--output"), help="choose the name for the output file WITHOUT extension"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-c", "--compare_by"), help="field of the metadata by which you want to compare")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

#q(save='no')


entropy = function(x)  {x = x[!is.na(x)]; x =x[x!=0]; p=x/sum(x); if (length(x) != 0) {return(-sum(p*log(p)))} else {NA}}

which_is_maj_iso = function(x) if (sum(x==0)==length(x)) {return(0)} else {return(which.max(x))}
gn_fun1 = function(x) max(table(as.numeric(x)))  # to be revised!!!
gn_fun2 = function(x) length(unique(x[x!=0]))


###########
## HUMAN ##
###########

# using only protein coding genes AND transcripts
human_gn_tx = read.table(opt$annotation, h=F, col.names=c('gene','tx'))
# transcript matrix
human_expr = read.table(opt$input_matrix, h=T)
human_expr[is.na(human_expr)] <- 0

# transcript matrix only for protein coding transcripts with gene id appended
human_tx_gn_expr = merge(human_gn_tx,human_expr,by.x='tx',by.y='row.names')
human_expr_iso = aggregate(human_tx_gn_expr[,-(1:2)], list(human_tx_gn_expr$gene), which_is_maj_iso)
rownames(human_expr_iso) = human_expr_iso[,1]
human_expr_iso[,1] = NULL
human_expr_iso_fun = data.frame(n_samples=apply(human_expr_iso,1,gn_fun1),n_maj_iso=apply(human_expr_iso,1,gn_fun2), organism='Human')

# add the number of annotated for human
hs_ann_iso = data.frame(table(human_gn_tx$gene))
names(hs_ann_iso) = c('gene','ann_iso')
human_expr_iso_fun_ann = merge(hs_ann_iso, human_expr_iso_fun, by.x='gene', by.y='row.names')

hs_mm_expr_iso_fun_ann = human_expr_iso_fun_ann


# read the metadata
mdata = read.table(opt$metadata, h=T, sep='\t')

variables = unique(mdata[,opt$compare_by])

shift = sapply(variables, function(x) {ids = unique(mdata[which(mdata[opt$compare_by]==x), 'labExpId']);
apply(human_expr_iso[ids], 1, function(x) if (length(unique(x))==1) {return(unique(x))}else{return(NA)} )
})

for (el in variables){
ids = unique(mdata[which(mdata[opt$compare_by]==el), 'labExpId'])
human_expr_iso[el] = apply(human_expr_iso[ids], 1, function(x) if (length(unique(x))==1) {return(unique(x))}else{return(NA)} )
}

shift2ggplot = data.frame(
# Total genes in the annotation
Total_genes = nrow(shift), 
# Genes with the same major isoform within the same group
same_maj_iso_within = nrow(shift[which(rowSums(shift==0) ==0),]),
# Genes with the same major isoform consistent within the group, but different between
diff_maj_iso_between = sum(apply(shift[which(rowSums(shift==0) ==0),], 1, function(x) length(unique(x)))>1) 
)

write.table(shift2ggplot, sprintf("%s.tsv", opt$output), sep='\t', quote=F, row.names=F)

## PLOT ##

library(ggplot2)
pdf(sprintf("%s.pdf", opt$output),width=12)

#ggplot(subset(hs_mm_expr_iso_fun_ann,ann_iso <= 25),aes(x=as.factor(ann_iso), y= n_samples)) + geom_boxplot() + facet_wrap(~organism,scale='free') + xlab('Number of annotated isoforms') + ylab('Max number of samples where the same major isoform is expressed') + opts(title='Most diffused major isoform')

ggplot(subset(hs_mm_expr_iso_fun_ann, ann_iso >= 2 & n_maj_iso>0),aes(x=as.factor(n_maj_iso),fill=cut(ann_iso, c(1,2,3,4,5,10,15,25,Inf)))) + geom_histogram(position='stack') + guides(fill = guide_legend(reverse = TRUE)) + scale_fill_hue(name='Total isoforms per gene',labels=c(2,3,4,5,'6-10','11-15','16-25','>25')) + labs(x='Major isoforms per gene', y='Number of protein coding genes') 

ggplot(subset(hs_mm_expr_iso_fun_ann,ann_iso<=35 & ann_iso>=2 & n_maj_iso>=1 ),aes(x=as.numeric(ann_iso))) + geom_histogram(aes(fill=cut(n_maj_iso, c(0,1,2,3,4,5,6,7,Inf))),position='stack',binwidth=1) + guides(fill=guide_legend(reverse=TRUE)) + stat_bin(aes(y=..count../x),breaks=seq(2,34,by=1),geom='line') + labs(x='Annotated isoforms per gene', y='Number of protein coding genes') + scale_fill_hue(name='Major isoforms per gene',labels=c(1,2,3,4,5,6,7,'>8'))

gp = ggplot(subset(hs_mm_expr_iso_fun_ann, ann_iso >= 2 & n_maj_iso>0), aes(x=ann_iso, y=n_maj_iso)) 
gp = gp + stat_bin2d(binwidth=c(1,1),aes(fill=cut(..count.., c(0,1,5,10,50,100,500,1000,2000,3000,Inf))),colour="black",size=.2) 
gp = gp + scale_fill_brewer(name='Number\nof genes', palette='Greens') + xlim(c(2,20)) 
gp = gp + ylim(c(0,20)) 
gp = gp + labs(x='Annotated isoforms per gene', y='Major isoforms per gene')
gp

#ggplot(subset(hs_mm_expr_iso_median_ann,hs_ann_iso!=1 & mm_ann_iso!=1 ),aes(x=floor(mm_entropy*10),y=floor(hs_entropy*10))) + stat_bin2d(binwidth=c(1,1),aes(fill=cut(..count.., c(0,1,5,10,100,500,1000,2000,3000,Inf))),colour ="black",size=.2) + scale_fill_brewer("count") + geom_abline(linetype=4) + opts(title="Median of entropy of isoform expression") + xlab('Mouse (entropy*10)') + ylab('Human (entropy*10)')
dev.off()
