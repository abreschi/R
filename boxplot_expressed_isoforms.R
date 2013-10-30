
## Q: How many isoforms are expressed per gene per cell line?

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-a", "--annotation"), help="two-column file with gene and tx ids"),
make_option(c("-o", "--output"), help="choose the name for the output file, WITHOUT extension"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-I", "--INPUT"), help="recover input if the table is already present", default=""),
make_option(c("-f", "--fill_by"), help="choose what to fill by [default=NA]", type='character', default=NA),
make_option(c("-c", "--color_by"), help="choose what to color by [default=NA]", type="character", default=NA)
)


parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

print(opt)

#--------------------
# FUNCTIONS
#--------------------

# FUNCTION 1
# This function calculates the relative expression of the most expressed isoform of a gene in a given condition.
# It returns NaN when the gene has no isoforms expressed
maj_iso_rel_expr = function(x) {max(x,na.rm=T)/sum(x,na.rm=T)}

# FUNCTION 2
# This function calculates the Shannon entropy for expressed isoforms of a gene in a given condition
# It returns NA when the gene has no isoforms expressed.
# The function is taken from the script:
source('~/R/functions.R')

# convert NAs to null
na2null = function(x) if(is.na(x)) {return(NULL)}else{return(x)}


###########
## HUMAN ##
###########

if (opt$INPUT == "") {
# (only whole cell PolyA+ is considered, ENCODE Jan11 two bioreplicates)
# gene, tx annotation file (only protein coding gene, but all txs)
human_gn_tx = read.table(opt$annotation, h=F, col.names=c('gene','tx'))
# transcript matrix
human_expr = read.table(opt$input_matrix, h=T)
# read the metadata
mdata = read.table(opt$metadata, h=T, sep="\t")
# replace NAs (bad IDR) with zeros
human_expr[is.na(human_expr)] <- 0

# add the gene locus to each tx
human_tx_gn_expr = merge(human_gn_tx, human_expr, by.x='tx', by.y='row.names')

#### EXPRESSED ISOFORMS
# count the number of expressed isoforms in each cell line
human_expr_iso = aggregate(human_tx_gn_expr[,-(1:2)], list(gene=human_tx_gn_expr$gene), function(x) sum(!is.na(x)&x!=0))
human_expr_iso_melt = melt(human_expr_iso, variable.name='sample_name', value.name='expr_iso')

#### MAJOR ISOFORM EXPRESSION
# calculate the relative epxression of the most expressed isoform for each gene in each sample
human_rel_expr_maj = aggregate(human_tx_gn_expr[,-(1:2)], list(gene=human_tx_gn_expr$gene), maj_iso_rel_expr)
human_rel_expr_maj_melt = melt(human_rel_expr_maj, variable.name='sample_name', value.name='rel_maj')

#### ISOFORM ENTROPY
# calculate the entropy of expressed isoforms for each gene in each sample
human_entr_iso = aggregate(human_tx_gn_expr[,-(1:2)], list(gene=human_tx_gn_expr$gene), entropy)
human_entr_iso_melt = melt(human_entr_iso, variable.name='sample_name', value.name='entr_iso')

#### ANNOTATED ISOFORMS
# count the number of annotated isoforms for each gene
human_ann_iso = as.data.frame(table(human_tx_gn_expr$gene))
names(human_ann_iso) = c('gene','ann_isoforms')

#### MERGING
# add the annotation information to the expressed isoforms
human_iso_all = merge(human_expr_iso_melt, human_rel_expr_maj_melt, by=c('gene', 'sample_name'))
human_iso_all = merge(human_iso_all, human_entr_iso_melt, by=c('gene', 'sample_name'))
human_iso_all = merge(human_iso_all, human_ann_iso, by = 'gene')
human_iso_all$organism = 'Human'

hs_mm_data_expr_ann_iso = human_iso_all
hs_mm_data_expr_ann_iso = merge(mdata[c("labExpId", opt$fill_by[!is.na(opt$fill_by)],  opt$color_by[!is.na(opt$color_by)])], 
hs_mm_data_expr_ann_iso, by.x='labExpId', by.y='sample_name')
write.table(hs_mm_data_expr_ann_iso, file=sprintf('%s.summary_isoform_expression.tsv',opt$output), sep='\t', quote=F, row.names=F)
}

if (opt$INPUT!="") {
hs_mm_data_expr_ann_iso = read.table(opt$INPUT, h=T)}

#q(save='no')
# label vectors with the number of genes having a given number of isoforms
#mm_lab = table(subset(hs_mm_data_expr_ann_iso, !duplicated(gene) & rnaExtract=='total')$ann_isoforms)
#head(subset(hs_mm_data_expr_ann_iso, !duplicated(gene) & tolower(rnaExtract)=='longpolya'))

hs_lab = table(subset(hs_mm_data_expr_ann_iso, !duplicated(gene) & organism == 'Human')$ann_isoforms)
#merged_labels = mapply(function(x,y) sprintf('%s\n%s', x,y), as.numeric(hs_lab)[1:30], as.numeric(mm_lab)[1:30])
merged_labels = as.numeric(hs_lab)[1:30]

##########
## plot ##
##########


# COMMENT: I have to add the number of genes in each boxplot

library(ggplot2)

theme_set(theme_bw(base_size=16))

pdf(sprintf("%s.pdf", opt$output), height=6, width=9)

max_iso = 25

gp = ggplot(subset(hs_mm_data_expr_ann_iso,ann_isoforms<=max_iso & expr_iso>0),aes(as.factor(ann_isoforms), expr_iso)) 
#gp = gp + geom_boxplot(fill='green') + ylim(c(0,20)) 
gp = gp + geom_boxplot(aes_string(fill=na2null(opt$fill_by), color=na2null(opt$color_by))) + ylim(c(0,max_iso)) 
gp = gp + labs(y='Number of isoforms expressed per gene', x='Number of annotated isoforms per gene') 
gp = gp + geom_abline(linetype='dotted', size=2, color='grey')
gp = gp + annotate('text', x = 1:max_iso, y = 0, label = names(hs_lab)[1:max_iso])
gp = gp + scale_fill_manual(values=cbbPalette)
gp = gp + scale_color_manual(values=rev(cbbPalette))
gp = gp + scale_x_discrete(labels=merged_labels) + theme(axis.text.x = element_text(angle=90, vjust=0.5))
#gp = gp + facet_wrap(~labExpId)
gp

gp = ggplot(subset(hs_mm_data_expr_ann_iso,ann_isoforms<=max_iso),aes(as.factor(ann_isoforms), rel_maj))
gp = gp + geom_boxplot(aes_string(fill=na2null(opt$fill_by), color=na2null(opt$color_by))) +ylim(c(0,1)) 
#gp = gp + geom_boxplot(fill='green') +ylim(c(0,1)) 
gp = gp + labs(y='Major isoform relative expression', x='Number of annotated isoforms per gene')
gp = gp + annotate('text', x = 1:max_iso, y = 0, label = names(hs_lab)[1:max_iso])
gp = gp + scale_fill_manual(values=cbbPalette)
gp = gp + scale_color_manual(values=rev(cbbPalette))
gp = gp + scale_x_discrete(labels=merged_labels) + theme(axis.text.x = element_text(angle=90, vjust=0.5))
gp

gp = ggplot(subset(hs_mm_data_expr_ann_iso,expr_iso<=15 & expr_iso>=2),aes(as.factor(expr_iso), entr_iso))
gp = gp + geom_boxplot(aes_string(fill=na2null(opt$fill_by), color=na2null(opt$color_by)))
gp = gp + labs(y = 'Shannon entropy', x = 'Number of expressed isoforms per gene')
gp = gp + geom_point(aes(x=as.factor(expr_iso), y = log(expr_iso+1)), color='red') 
#gp = gp + annotate('text', x = (2:max_iso)-1, y = -1/4, label = names(hs_lab)[2:max_iso])
gp = gp + scale_fill_manual(values=cbbPalette)
gp = gp + scale_color_manual(values=rev(cbbPalette))
gp

gp = ggplot(subset(hs_mm_data_expr_ann_iso,ann_isoforms<=max_iso & ann_isoforms>=2),aes(as.factor(ann_isoforms), entr_iso))
gp = gp + geom_boxplot(aes_string(fill=na2null(opt$fill_by), color=na2null(opt$color_by)))
#gp = gp + geom_boxplot(fill='green') 
gp = gp + labs(y = 'Shannon entropy', x = 'Number of annotated isoforms per gene')
gp = gp + geom_point(aes(x=as.factor(ann_isoforms), y = log(ann_isoforms+1)), color='red') 
gp = gp + annotate('text', x = (2:max_iso)-1, y = -1/4, label = names(hs_lab)[2:max_iso])
gp = gp + scale_fill_manual(values=cbbPalette)
gp = gp + scale_color_manual(values=rev(cbbPalette))
gp

dev.off()

ggsave(sprintf("%s.jpeg", opt$output), h=6, w=9)

