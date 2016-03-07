#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


##################
# OPTION PARSING
##################


option_list <- list(

make_option(c("-i", "--input_matrix"), 
	help="the matrix you want to analyze"),

make_option(c("-a", "--annotation"), 
	help="two-column file with gene and tx ids, no header."),

make_option(c("-o", "--output"),
	help="choose the name for the output file, WITHOUT extension"),

make_option(c("-m", "--metadata"),
	help="tsv file with metadata on matrix experiment"),

make_option(c("-f", "--fill_by"),
	help="choose what to fill by. Leave empty for no filling"),

make_option(c("-c", "--color_by"),
	help="choose what to color by. Leave empty for no coloring")

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
source('~abreschi/R/functions.R')



###########
## HUMAN ##
###########


# gene, tx annotation file (only protein coding gene, but all txs)
human_gn_tx = read.table(opt$annotation, h=F, col.names=c('gene','tx'))

# transcript matrix
human_expr = read.table(opt$input_matrix, h=T)

# read the metadata
mdata = read.table(opt$metadata, h=T, sep="\t")
mdata[,"labExpId"] <- gsub(",", ".", mdata[,"labExpId"]) 

# replace NAs (bad IDR) with zeros
human_expr[is.na(human_expr)] <- 0

# add the gene locus to each tx
human_tx_gn_expr = merge(human_gn_tx, human_expr, by.x='tx', by.y='row.names')


#### EXPRESSED ISOFORMS
# count the number of expressed isoforms in each cell line
human_expr_iso = aggregate(human_tx_gn_expr[,-(1:2)], list(gene=human_tx_gn_expr$gene), function(x) sum(!is.na(x)&x!=0))
# equivalent to:
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
human_ann_iso = setNames(aggregate(tx~gene, human_tx_gn_expr, length), c("gene", "ann_isoforms"))

#### MERGING
# add the annotation information to the expressed isoforms
human_iso_all = merge(human_expr_iso_melt, human_rel_expr_maj_melt, by=c('gene', 'sample_name'))
human_iso_all = merge(human_iso_all, human_entr_iso_melt, by=c('gene', 'sample_name'))
human_iso_all = merge(human_iso_all, human_ann_iso, by = 'gene')

hs_mm_data_expr_ann_iso = human_iso_all
hs_mm_data_expr_ann_iso = merge(unique(mdata[c("labExpId", opt$fill_by, opt$color_by)]),
hs_mm_data_expr_ann_iso, by.x='labExpId', by.y='sample_name')
write.table(hs_mm_data_expr_ann_iso, file=sprintf('%s.summary_isoform_expression.tsv',opt$output), sep='\t', quote=F, row.names=F)


merged_labels = aggregate(gene~tx, aggregate(tx~gene, human_tx_gn_expr, length), length)

##########
## plot ##
##########


# COMMENT: I have to add the number of genes in each boxplot

theme_set(theme_bw(base_size=16))

pdf(sprintf("%s.pdf", opt$output), height=6, width=9)

max_iso = 15

data = subset(hs_mm_data_expr_ann_iso,ann_isoforms<=max_iso & expr_iso>0)
gp = ggplot(data, aes(as.factor(ann_isoforms), expr_iso)) 
gp = gp + geom_boxplot(aes_string(fill=opt$fill_by, color=opt$color_by)) 
gp = gp + ylim(c(0, max_iso)) 
gp = gp + labs(y='Number of isoforms expressed per gene', x='Number of annotated isoforms per gene') 
gp = gp + geom_abline(linetype='dotted', size=2, color='grey')
gp = gp + geom_text(data=subset(merged_labels, tx<=max_iso), aes(x = tx, y = 0, label = gene), angle=90)
gp = gp + scale_fill_manual(values=cbbPalette)
gp

data = subset(hs_mm_data_expr_ann_iso,ann_isoforms<=max_iso)
gp = ggplot(data, aes(as.factor(ann_isoforms), rel_maj))
gp = gp + geom_boxplot(aes_string(fill=opt$fill_by, color=opt$color_by))
gp = gp + ylim(c(0,1)) 
gp = gp + labs(y='Major isoform relative expression', x='Number of annotated isoforms per gene')
gp = gp + geom_text(data=subset(merged_labels, tx<=max_iso), aes(x = tx, y = 0, label = gene), angle=90)
gp = gp + scale_fill_manual(values=cbbPalette)
gp = gp + scale_x_discrete(labels=merged_labels) + theme(axis.text.x = element_text(angle=90, vjust=0.5))
gp

data = subset(hs_mm_data_expr_ann_iso,expr_iso<=15 & expr_iso>=2)
gp = ggplot(data, aes(as.factor(expr_iso), entr_iso))
gp = gp + geom_boxplot(aes_string(fill=opt$fill_by, color=opt$color_by))
gp = gp + labs(y = 'Shannon entropy', x = 'Number of expressed isoforms per gene')
gp = gp + geom_point(aes(x=as.factor(expr_iso), y = log(expr_iso+1)), color='red') 
gp = gp + scale_fill_manual(values=cbbPalette)
gp

data = subset(hs_mm_data_expr_ann_iso,ann_isoforms<=max_iso & ann_isoforms>=2)
gp = ggplot(data, aes(as.factor(ann_isoforms), entr_iso))
gp = gp + geom_boxplot(aes_string(fill=opt$fill_by, color=opt$color_by))
#gp = gp + ylim(c(0, max_iso)) 
gp = gp + labs(y = 'Shannon entropy', x = 'Number of annotated isoforms per gene')
gp = gp + geom_point(aes(x=as.factor(ann_isoforms), y = log(ann_isoforms+1)), color='red') 
gp = gp + geom_text(data=subset(merged_labels, tx<=max_iso), aes(x = tx, y = 0, label = gene), angle=90)
gp = gp + scale_fill_manual(values=cbbPalette)
gp

dev.off()

ggsave(sprintf("%s.jpeg", opt$output), h=6, w=9)

