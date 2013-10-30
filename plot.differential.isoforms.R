
## Q: How many isoforms are expressed per gene per cell line?

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix with the trancsript rpkms"),
make_option(c("-a", "--annotation"), help="tsv with gene_id, tx_id, gene_type, gene_name, transcript_name IN THIS ORDER"),
make_option(c("-o", "--output"), help="choose the name for the output file"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-d", "--de"), help='output of edgeR for differentially expressed genes')

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

# convert NAs to nul
na2null = function(x) if(is.na(x)) {return(NULL)}else{return(x)}


###########
## HUMAN ##
###########

# (only whole cell PolyA+ is considered, ENCODE Jan11 two bioreplicates)
# gene, tx annotation file (only protein coding gene, but all txs)
human_gn_tx = read.table(opt$annotation, h=F, col.names=c('gene','tx','gene_type','gene_name','transcript_type'))
human_gn_tx = (merge(human_gn_tx, as.data.frame(table(human_gn_tx$gene), responseName='ann_iso'), by.x='gene', by.y='Var1'))

# transcript matrix
human_expr = read.table(opt$input_matrix, h=T)
# read the metadata
mdata = read.table(opt$metadata, h=T)
# replace NAs (bad IDR) with zeros
human_expr[is.na(human_expr)] <- 0

# add the gene locus to each tx
human_tx_gn_expr = merge(human_gn_tx, human_expr, by.x='tx', by.y='row.names')

# DE table
DE = read.table(opt$de, h=T)
DE = subset(merge(human_gn_tx, DE, by.y='row.names', by.x='tx'), FDR<0.01)
DE = merge(DE, as.data.frame(table(DE$gene), responseName='de_tx'), by.y='Var1', by.x='gene')

ggplot(subset(DE, ann_iso<=25)) + geom_boxplot(aes(x=factor(ann_iso), y=de_tx) )


# shown in the slides
case = subset(human_tx_gn_expr, gene=='ENSG00000035681')

df = merge(melt(case), unique(mdata[c('labExpId','cell')]), by.x = 'variable', by.y='labExpId')
df = ddply(df, .(tx,cell), transform, avg=mean(value), sdev = sd(value))

gp = ggplot(subset(df, tx %in% c('ENST00000038176', 'ENST00000523177', 'ENST00000519174', 'ENST00000523106', 'ENST00000521972')), 
aes(x=cell, y=avg))
gp = gp + geom_line(aes(group=tx, color=tx, linetype=transcript_type)) 
gp + geom_errorbar(aes(ymin=avg-sdev, ymax=avg+sdev, color=tx), width=.02)


# for testing
case = subset(human_tx_gn_expr, gene=='ENSG00000010017')

df = merge(melt(case), unique(mdata[c('labExpId','cell')]), by.x = 'variable', by.y='labExpId')
df = ddply(df, .(tx,cell), transform, avg=mean(value), sdev = sd(value))

gp = ggplot(df, aes(x=cell, y=avg))
gp = gp + geom_line(aes(group=tx, color=tx, linetype=transcript_type)) 
gp + geom_errorbar(aes(ymin=avg-sdev, ymax=avg+sdev, color=tx), width=.02)



case = subset(human_tx_gn_expr, gene=='ENSG00000162517')

df = merge(melt(case), unique(mdata[c('labExpId','cell')]), by.x = 'variable', by.y='labExpId')
df = ddply(df, .(tx,cell), transform, avg=mean(value), sdev = sd(value))

gp = ggplot(df,aes(x=cell, y=avg))
gp = gp + geom_line(aes(group=tx, color=tx, linetype=transcript_type)) 
gp + geom_errorbar(aes(ymin=avg-sdev, ymax=avg+sdev, color=tx), width=.02)


DE_no_genes = subset(DE, !(gene %in% genes$merged)) 


genes = head(subset(DE, FDR<0.01 & ann_iso==4 & abs(logFC)>3)$gene,3)

df = merge(melt(subset(human_tx_gn_expr, gene %in% genes)), unique(mdata[c('labExpId','cell')]), by.x = 'variable', by.y='labExpId')
df = ddply(df, .(tx,cell), transform, avg=mean(value), sdev=sd(value))


gp = ggplot(df, aes(x=cell, y=avg)) 
gp = gp + geom_line(aes(group=tx, color=tx, linetype=transcript_type)) 
gp = gp + facet_wrap(~gene_name, scale='free_y', nrow=1)  
gp + geom_errorbar(aes(ymin=avg-sdev, ymax=avg+sdev, color=tx), width=0.07) 









q(save='no')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

pdf(opt$output, height=7, width=10)

gp = ggplot(subset(hs_mm_data_expr_ann_iso,ann_isoforms<=20 & expr_iso>0),aes(as.factor(ann_isoforms), expr_iso)) 
#gp = gp + geom_boxplot(fill='green') + ylim(c(0,20)) 
gp = gp + geom_boxplot(aes_string(fill=na2null(opt$fill_by), color=na2null(opt$color_by))) + ylim(c(0,20)) 
gp = gp + labs(y='Number of isoforms expressed per gene', x='Number of annotated isoforms per gene') 
gp = gp + geom_abline(linetype='dotted', size=2, color='grey')
gp = gp + annotate('text', x = 1:20, y = 0, label = names(hs_lab)[1:20])
gp = gp + scale_fill_manual(values=cbbPalette)
gp = gp + scale_color_manual(values=rev(cbbPalette))
gp = gp + scale_x_discrete(labels=merged_labels) + theme(axis.text.x = element_text(angle=90, vjust=0.5))
#gp = gp + facet_wrap(~labExpId)
gp

gp = ggplot(subset(hs_mm_data_expr_ann_iso,ann_isoforms<=20),aes(as.factor(ann_isoforms), rel_maj))
gp = gp + geom_boxplot(aes_string(fill=na2null(opt$fill_by), color=na2null(opt$color_by))) +ylim(c(0,1)) 
#gp = gp + geom_boxplot(fill='green') +ylim(c(0,1)) 
gp = gp + labs(y='Major isoform relative expression', x='Number of annotated isoforms per gene')
gp = gp + annotate('text', x = 1:20, y = 0, label = names(hs_lab)[1:20])
gp = gp + scale_fill_manual(values=cbbPalette)
gp = gp + scale_color_manual(values=rev(cbbPalette))
gp = gp + scale_x_discrete(labels=merged_labels) + theme(axis.text.x = element_text(angle=90, vjust=0.5))
gp

gp = ggplot(subset(hs_mm_data_expr_ann_iso,ann_isoforms<=20 & ann_isoforms>=2),aes(as.factor(ann_isoforms), entr_iso))
gp = gp + geom_boxplot(aes_string(fill=na2null(opt$fill_by), color=na2null(opt$color_by)))
#gp = gp + geom_boxplot(fill='green') 
gp = gp + labs(y = 'Shannon entropy', x = 'Number of annotated isoforms per gene')
gp = gp + geom_point(aes(x=as.factor(ann_isoforms), y = log(ann_isoforms+1)), color='red') 
gp = gp + annotate('text', x = (2:20)-1, y = -1/4, label = names(hs_lab)[2:20])
gp = gp + scale_fill_manual(values=cbbPalette)
gp = gp + scale_color_manual(values=rev(cbbPalette))
gp

dev.off()
