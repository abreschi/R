# This scripts identifies tissue-specific genes based on logratio between the two highest rpkms


##------------
## LIBRARIES
##------------ 

cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library("optparse"))
cat("DONE\n\n")


options(stringsAsFactor=FALSE)
theme_set(theme_bw())


##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), 
	help="the matrix you want to analyze"),

make_option(c("-m", "--metadata"), 
	help="tsv file with metadata on matrix experiment. If empty, column names are used as labels."),

make_option(c("-o", "--output"), default="out",
	help="additional tags for otuput [default=%default]"),

make_option(c("-b", "--breaks"), type="double", default=3, 
	help="threshold on the breaks for plotting [default=%default]"),

make_option(c("-p", "--pval"), type="double", default=0.01, 
	help="p-value threshold for plotting [default=%default]"),

make_option(c("-t", "--tags"), default="labExpId",
	help="comma-separated field names you want to display in the labels. Leave default for using column names [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="verbose output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

##############
# BEGIN
##############

m = read.table(opt$input_matrix, h=T)

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}



###############################
## LOG-FOLD MULTIPLE SAMPLES ##
###############################

#---------------
# Functions
#--------------

# calculate the max distance between two elements in a list
max_diff = function(x)  {
	if (sum(x,na.rm=T) == 0) {                                   # if the array is all 0s or mixed 0s and NAs, return NA (no values)
		return(data.frame("log10_max_diff"=NA_real_, "break_max_diff"=NA_integer_, "rgt_max_diff"="",stringsAsFactors=F))}
	if (sum(!is.na(x) & x!= 0) == 1) {                                   # if there is only one value, maxdiff/sd is the sqrt(N)
		return(data.frame("log10_max_diff"=sqrt(length(x)), "break_max_diff"=length(x)-1, "rgt_max_diff"=names(which.max(x)),stringsAsFactors=F))}
	x[x == 0] <- NA
    x_sorted = sort(log10(x), na.last=F)                              # sort the array, apply log10 and remove 0s
	# if the array has all NAs but one, diff() returns an array with all NAs, 
    # then which.max() will return an array long 0, that makes indexing crash
	diffs = diff(as.numeric(x_sorted))
	log10_max_diff = max(diffs, na.rm=T)/sd(as.numeric(x_sorted),na.rm=T)       # find the maximum difference
	break_max_diff = which.max(diffs)
	rgt_max_diff = names(x_sorted)[-(1:break_max_diff)]
	return(data.frame("log10_max_diff"=log10_max_diff, "break_max_diff"=break_max_diff, "rgt_max_diff"=rgt_max_diff, stringsAsFactors=F))
	}



# Find the 1-cdf from an estimated density function
pval = function(x, pdf) {if (is.na(x))return(NA); f=approxfun(pdf$x, pdf$y, yleft=0, yright=0); p=integrate(f,x,Inf,rel.tol=0.001)$value; 
	return(p)}

# shuffle all cells in the matrix
mat_shuffle = function(m) {set.seed(123)
df<-data.frame(matrix(sample(unlist(m), length(unlist(m))), nrow=nrow(m), ncol=ncol(m))); 
colnames(df)<-colnames(m); rownames(df)<-rownames(m); return(df)}


# I want to evaluate the impact of zeros and NAs since I have to do the log
# Problems:
# - all NAs but one --> no stdev. I assume all the other values are equal and I assign sqrt(N)
# - at least one zero --> it will alter all the other diff because I would get an Inf ratio, so I treat it as NA. 
#                         The sd will have less points
# - all zeros but one --> Inf, so I treat all the zeros as NAs, and I assing sqrt(N)
# - all zeros or NAs --> return NA for everything


#----------#
# BEGIN    #  
#----------#

maxdiff <- do.call(rbind, apply(m, 1, max_diff))
maxdiff$gene <- sapply(rownames(maxdiff), function(x) strsplit(x, ".",fixed=T)[[1]][1])

matshuff <- mat_shuffle(m)

maxdiff_matshuff <- do.call(rbind, apply(matshuff, 1, max_diff))
maxdiff_matshuff$gene <- sapply(rownames(maxdiff_matshuff), function(x) strsplit(x, ".",fixed=T)[[1]][1])

maxdiff_matshuff_pdf <- density(subset(maxdiff_matshuff, !duplicated(gene))$log10_max_diff,na.rm=T)

maxdiff$pval <- sapply(maxdiff$log10_max_diff, pval, pdf=maxdiff_matshuff_pdf)
#maxdiff$pvcor <- sapply(maxdiff$pval, function(x) 1-((1-x)**nrow(m)))
maxdiff$pval.adj = p.adjust(maxdiff$pval, method="BH")

maxdiff_matshuff$pval <- sapply(maxdiff_matshuff$log10_max_diff, pval, pdf=maxdiff_matshuff_pdf)
#maxdiff_matshuff$pvcor <- sapply(maxdiff_matshuff$pval, function(x) 1-((1-x)**nrow(m)))
maxdiff_matshuff$pval.adj = p.adjust(maxdiff_matshuff$pval, method="BH")

thr = quantile(subset(maxdiff_matshuff, !duplicated(gene))$log10_max_diff, 0.95)


mdata = read.table(opt$metadata, h=T, sep='\t')
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))

tags = strsplit(opt$tags, ",")[[1]]
mdata = unique(mdata[unique(c("labExpId",tags))])
maxdiff <- merge(maxdiff, mdata, by.x="rgt_max_diff", by.y="labExpId") 
maxdiff_matshuff <- merge(maxdiff_matshuff, unique(mdata[unique(c("labExpId",tags))]), by.x="rgt_max_diff", by.y="labExpId") 
maxdiff$label = apply(maxdiff[tags], 1, paste, collapse=".")
maxdiff_matshuff$label = apply(maxdiff_matshuff[tags], 1, paste, collapse=".")

# ****

# Print output
output = sprintf("maxdiff.%s", opt$output)
write.table(maxdiff, sprintf("%s.tsv", output), row.names=F, quote=F, sep="\t")
write.table(maxdiff_matshuff, sprintf("%s.rand.tsv", output), row.names=F, quote=F, sep="\t")

# ****

# Diagnostic plots
maxdiff$status = "obs"; maxdiff_matshuff$status = "pred"
df = subset(rbind(maxdiff,maxdiff_matshuff), !duplicated(interaction(gene,status)))

pdf(sprintf("%s.pdf", output), w=10)

gp = ggplot(df, aes(x=log10_max_diff)) 
gp = gp + geom_histogram(aes(fill=status), color='black', position='dodge')
gp = gp + scale_fill_manual(values=c("brown", "orange"))
gp

gp = ggplot(df, aes(x=as.factor(break_max_diff)))
gp = gp + geom_histogram(aes(fill=status), color='black', position='dodge')
gp = gp + scale_fill_manual(values=c("blue", "cyan"))
gp

gp = ggplot(df, aes(x=as.factor(break_max_diff)))
gp = gp + geom_boxplot(aes(y=log10_max_diff, fill=status))
gp = gp + scale_fill_manual(values=c("darkgreen", "green"))
gp

maxdiff = ddply(maxdiff, .(label), transform, overcounts=sum(pval<opt$pval & break_max_diff>=ncol(m)-opt$breaks))
maxdiff$label <- factor(maxdiff$label, levels=unique(maxdiff[order(maxdiff$overcounts, decreasing=T),"label"]))

maxdiff_matshuff = ddply(maxdiff_matshuff, .(label), transform, overcounts=sum(pval<opt$pval & break_max_diff>=ncol(m)-opt$breaks))
maxdiff_matshuff$label = factor(maxdiff_matshuff$label, levels=unique(maxdiff[order(maxdiff$overcounts, decreasing=T),"label"]))

gp = ggplot(ddply(maxdiff, .(label), summarise, counts=sum(pval<opt$pval & break_max_diff>=ncol(m)-opt$breaks)), aes(x=label)) 
gp = gp + geom_bar(color='black',aes(fill=counts,y=counts), stat="identity")
gp = gp + scale_fill_gradientn(colours = rev(heat.colors(16)))
gp = gp + geom_bar(data=ddply(maxdiff_matshuff, .(label), summarise, counts=sum(pval<opt$pval & break_max_diff>=ncol(m)-opt$breaks)),
fill='black',alpha=0.5, aes(x=label,y=counts), stat="identity")
gp = gp + theme(axis.text.x = element_text(angle=45, hjust=1))
gp = gp + labs(title=sprintf("Number of tissue_specific over-expressed genes\nbreak_max_diff>=%s, pval<=%s",ncol(m)-opt$breaks,opt$pval))
gp


maxdiff = ddply(maxdiff, .(label), transform, undercounts=sum(pval<opt$pval & break_max_diff<=opt$breaks))
maxdiff$label <- factor(maxdiff$label, levels=unique(maxdiff[order(maxdiff$undercounts, decreasing=T),"label"]))

maxdiff_matshuff = ddply(maxdiff_matshuff, .(label), transform, undercounts=sum(pval<opt$pval & break_max_diff<=opt$breaks))
maxdiff_matshuff$label = factor(maxdiff_matshuff$label, levels=unique(maxdiff[order(maxdiff$undercounts, decreasing=T),"label"]))

gp = ggplot(ddply(maxdiff, .(label), summarise, counts=sum(pval<opt$pval & break_max_diff<=opt$breaks)), aes(x=label)) 
gp = gp + geom_bar(color='black',aes(fill=counts,y=counts), stat="identity")
gp = gp + scale_fill_gradientn(colours = rev(heat.colors(16)))
gp = gp + geom_bar(data=ddply(maxdiff_matshuff, .(label), summarise, counts=sum(pval<opt$pval & break_max_diff<=opt$breaks)),
fill='black',alpha=0.5, aes(x=label,y=counts), stat="identity")
gp = gp + theme(axis.text.x = element_text(angle=45, hjust=1))
gp = gp + labs(title=sprintf("Number of tissue_specific under-expressed genes\nbreak_max_diff<=%s, pval<=%s", opt$breaks,opt$pval))
gp

dev.off()

q(save='no')








#%%% hsCSHL %%%#
hsCSHL_maxdiff <- ddply(hsCSHL_orth, c('hs','mm'), function(x) max_diff(x[-(1:2)]))

# generate the random sets or human data
hsCSHL_rnorm <- cbind(orth, rnorm_df(nrow(hsCSHL_orth), ncol(hsCSHL_orth)-2, lab=names(hsCSHL_orth[-(1:2)])))
hsCSHL_colshuf <- col_shuffle(hsCSHL_orth)
hsCSHL_matshuf <- mat_shuffle(hsCSHL_orth, c('hs','mm')) 

# do the computations for the random samples
hsCSHL_rnorm_maxdiff <- ddply(hsCSHL_rnorm, c('hs','mm'), function(x) max_diff_rnorm(x[-(1:2)]))
hsCSHL_colshuf_maxdiff <- ddply(hsCSHL_colshuf, c('hs','mm'), function(x) max_diff(x[-(1:2)]))
hsCSHL_matshuf_maxdiff <- ddply(hsCSHL_matshuf, c('hs','mm'), function(x) max_diff(x[-(1:2)]))

hsCSHL_rnorm_maxdiff_pdf <- density(subset(hsCSHL_rnorm_maxdiff, !duplicated(hs))$log10_max_diff,na.rm=T)
hsCSHL_colshuf_maxdiff_pdf <- density(subset(hsCSHL_colshuf_maxdiff, !duplicated(hs))$log10_max_diff,na.rm=T)
hsCSHL_matshuf_maxdiff_pdf <- density(subset(hsCSHL_matshuf_maxdiff, !duplicated(hs))$log10_max_diff,na.rm=T)

hsCSHL_maxdiff$pval <- sapply(hsCSHL_maxdiff$log10_max_diff, pval, pdf = hsCSHL_matshuf_maxdiff_pdf)
hsCSHL_maxdiff$pvcor <- sapply(hsCSHL_maxdiff$pval, function(x) 1-((1-x)**nrow(orth)))

hsCSHL_all = rbind.fill(cbind(hsCSHL_maxdiff, dataset='hsCSHL', sampling='real'),
	cbind(hsCSHL_rnorm_maxdiff, dataset='hsCSHL', sampling='rnorm'),
	cbind(hsCSHL_colshuf_maxdiff, dataset='hsCSHL', sampling='colshuf'),
	cbind(hsCSHL_matshuf_maxdiff, dataset='hsCSHL', sampling='matshuf'))

hsCSHL_thr = quantile(subset(hsCSHL_rnorm_maxdiff, !duplicated(hs))$log10_max_diff, 0.95)
#hsCSHL_thr_bonf = quantile(subset(hsCSHL_rnorm_maxdiff, !duplicated(hs))$log10_max_diff, (0.95)**(1/nrow(orth))  )

ggplot(subset(hsCSHL_all, !duplicated(paste(hsCSHL_all$hs,hsCSHL_all$dataset, hsCSHL_all$sampling))), aes(x=log10_max_diff)) + 
    geom_density(aes(col=sampling)) + geom_vline(x=hsCSHL_thr) +
    geom_vline(x=hsCSHL_thr_bonf)
ggplot(subset(hsCSHL_all, !duplicated(paste(hsCSHL_all$hs,hsCSHL_all$dataset, hsCSHL_all$sampling))), aes(x=break_max_diff)) + 
    geom_density(aes(col=sampling))




#%%% mmCSHL %%%#

mm_maxdiff <- ddply(mmCSHL_orth, c('hs','mm'), function(x) max_diff(x[-(1:2)]))

# generate the random sets or human data
mm_rnorm <- cbind(orth, rnorm_df(nrow(mmCSHL_orth), ncol(mmCSHL_orth)-2, lab=names(mmCSHL_orth[-(1:2)])))
mm_colshuf <- col_shuffle(mmCSHL_orth)
mm_matshuf <- mat_shuffle(mmCSHL_orth, c('hs','mm')) 

# do the computations for the random samples
mm_rnorm_maxdiff <- ddply(mm_rnorm, c('hs','mm'), function(x) max_diff_rnorm(x[-(1:2)]))
mm_colshuf_maxdiff <- ddply(mm_colshuf, c('hs','mm'), function(x) max_diff(x[-(1:2)]))
mm_matshuf_maxdiff <- ddply(mm_matshuf, c('hs','mm'), function(x) max_diff(x[-(1:2)]))

mm_rnorm_maxdiff_pdf <- density(subset(mm_rnorm_maxdiff, !duplicated(hs))$log10_max_diff,na.rm=T)
mm_colshuf_maxdiff_pdf <- density(subset(mm_colshuf_maxdiff, !duplicated(hs))$log10_max_diff,na.rm=T)
mm_matshuf_maxdiff_pdf <- density(subset(mm_matshuf_maxdiff, !duplicated(hs))$log10_max_diff,na.rm=T)

mm_maxdiff$pval <- sapply(mm_maxdiff$log10_max_diff, pval, pdf = mm_matshuf_maxdiff_pdf)
mm_maxdiff$pvcor <- sapply(mm_maxdiff$pval, function(x) 1-((1-x)**nrow(orth)))

mm_thr = quantile(subset(mm_rnorm_maxdiff, !duplicated(hs))$log10_max_diff, 0.95)

mm_all = rbind.fill(cbind(mm_maxdiff, dataset='mm', sampling='real'),
	cbind(mm_rnorm_maxdiff, dataset='mm', sampling='rnorm'),
	cbind(mm_colshuf_maxdiff, dataset='mm', sampling='colshuf'),
	cbind(mm_matshuf_maxdiff, dataset='mm', sampling='matshuf'))


#%%% hsHBM %%%#

hsHBM_maxdiff <- ddply(hsHBM_orth, c('hs','mm'), function(x) max_diff(x[-(1:2)]))

# generate the random sets or human data
hsHBM_rnorm <- cbind(orth, rnorm_df(nrow(hsHBM_orth), ncol(hsHBM_orth)-2, lab=names(hsHBM_orth[-(1:2)])))
hsHBM_colshuf <- col_shuffle(hsHBM_orth)
hsHBM_matshuf <- mat_shuffle(hsHBM_orth, c('hs','mm')) 

# do the computations for the random samples
hsHBM_rnorm_maxdiff <- ddply(hsHBM_rnorm, c('hs','mm'), function(x) max_diff_rnorm(x[-(1:2)]))
hsHBM_colshuf_maxdiff <- ddply(hsHBM_colshuf, c('hs','mm'), function(x) max_diff(x[-(1:2)]))
hsHBM_matshuf_maxdiff <- ddply(hsHBM_matshuf, c('hs','mm'), function(x) max_diff(x[-(1:2)]))

hsHBM_rnorm_maxdiff_pdf <- density(subset(hsHBM_rnorm_maxdiff, !duplicated(hs))$log10_max_diff,na.rm=T)
hsHBM_colshuf_maxdiff_pdf <- density(subset(hsHBM_colshuf_maxdiff, !duplicated(hs))$log10_max_diff,na.rm=T)
hsHBM_matshuf_maxdiff_pdf <- density(subset(hsHBM_matshuf_maxdiff, !duplicated(hs))$log10_max_diff,na.rm=T)

hsHBM_maxdiff$pval <- sapply(hsHBM_maxdiff$log10_max_diff, pval, pdf = hsHBM_matshuf_maxdiff_pdf)
hsHBM_maxdiff$pvcor <- sapply(hsHBM_maxdiff$pval, function(x) 1-((1-x)**nrow(orth)))

hsHBM_thr = quantile(subset(hsHBM_rnorm_maxdiff, !duplicated(hs))$log10_max_diff, 0.95)

hsHBM_all = rbind.fill(cbind(hsHBM_maxdiff, dataset='hsHBM', sampling='real'),
	cbind(hsHBM_rnorm_maxdiff, dataset='hsHBM', sampling='rnorm'),
	cbind(hsHBM_colshuf_maxdiff, dataset='hsHBM', sampling='colshuf'),
	cbind(hsHBM_matshuf_maxdiff, dataset='hsHBM', sampling='matshuf'))


# Find common genes for the clustering

ALL = rbind.fill(hsCSHL_all, mm_all, hsHBM_all)
#ALL_thr = rbind.fill(subset(hsCSHL_maxdiff, log10_max_diff >= hsCSHL_thr & break_max_diff>=30), 
#	subset(mm_maxdiff, log10_max_diff >= mm_thr &break_max_diff>=25),
#	subset(hsHBM_maxdiff, log10_max_diff >= hsHBM_thr &break_max_diff >=10))

ALL_thr = rbind.fill(subset(hsCSHL_all, log10_max_diff >= hsCSHL_thr & break_max_diff>=30), 
	subset(mm_all, log10_max_diff >= mm_thr &break_max_diff>=25),
	subset(hsHBM_all, log10_max_diff >= hsHBM_thr &break_max_diff >=10))
# to Dong
write.table(subset(ALL_thr, dataset=='mm' & sampling=='real'&rgt_max_diff=='mouse_Testis_adult.8wks_CSHL')[,c('hs','mm')],
 quote=F,file='mm_testis.tsv',sep='\t' , row.names=F, col.names=F)
# write the tissue-specific genes
write.table(ALL_thr, quote=F, file='maxdiff_tissue-specific.tsv', sep='\t', row.names=F, col.names=F)

#ALL_thr = rbind.fill(subset(hsCSHL_all, log10_max_diff >= hsCSHL_thr), 
#	subset(mm_all, log10_max_diff >= mm_thr),
#	subset(hsHBM_all, log10_max_diff >= hsHBM_thr))

pdf('maxdiff_distribution.pdf', width=15, height=10)
gp = ggplot(subset(ALL, !duplicated(paste(ALL$hs, ALL$dataset, ALL$sampling))), aes(x=log10_max_diff)) 
gp = gp + geom_density(aes(color=sampling)) + scale_color_manual(values=cbbPalette)
gp = gp + facet_wrap(~dataset, scale='free', ncol=2)
gp = gp + theme(legend.position=c(0.6, 0.4)) + labs(x='max_diff_log10 / sd')
gp #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gp = ggplot(subset(ALL, !duplicated(paste(ALL$hs, ALL$dataset, ALL$sampling))), aes(x=break_max_diff)) 
gp = gp + geom_density(aes(color=sampling)) + scale_color_manual(values=cbbPalette)
gp = gp + facet_wrap(~dataset, scale='free', ncol=2)
gp = gp + theme(legend.position=c(0.6, 0.4))
gp #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gp = ggplot(ALL_thr , aes(x=rgt_max_diff)) 
gp = gp + geom_histogram(aes(fill=sampling), alpha=1/4, position = 'identity')
gp = gp + facet_wrap(~dataset, scale='free', ncol=2) + scale_fill_manual(values=cbbPalette)
gp = gp + theme(axis.text.x = element_text(angle=90), legend.position=c(0.6, 0)) 
gp = gp + labs(x='', y='Number of genes')
gp
dev.off()



ALL_real_thr = rbind.fill(subset(hsCSHL_maxdiff, log10_max_diff >= hsCSHL_thr & break_max_diff>=30), 
	subset(mm_maxdiff, log10_max_diff >= mm_thr &break_max_diff>=25),
	subset(hsHBM_maxdiff, log10_max_diff >= hsHBM_thr &break_max_diff >=10))



common_molten = expand.grid(x= names(hsCSHL_hsHBM_mmCSHL_orth)[-(1:2)], y=names(hsCSHL_hsHBM_mmCSHL_orth)[-(1:2)], stringsAsFactors=F)
common_molten$common = mapply( function(a,b) {set1=as.character(subset(ALL_real_thr, rgt_max_diff==a)$hs);
	set2=as.character(subset(ALL_real_thr, rgt_max_diff==b)$hs);
#	return(length(intersect(set1,set2)))},
	return(length(intersect(set1,set2))/min(length(set1), length(set2)))},
	common_molten$x, common_molten$y)
common_unmelt = as.matrix(cast(common_molten, x~y))

ALL_real_thr_expr = merge(unique(ALL_real_thr[,c('hs','mm')]), hsCSHL_hsHBM_mmCSHL_orth, by=c('hs','mm'))


pdf('maxdiff_cluster_tissueSpecific.pdf', height=30, w=22)
#side colors for rows
my_colors = list('human_CSHL'='springgreen1', 'human_HBM'='steelblue1', 'mouse'='violet')
sample_names = rownames(common_unmelt)
rc = rep(0,length(sample_names))
rc[grep('mouse', sample_names)] = as.character(my_colors['mouse'])
rc[grep('HBM', sample_names)] = as.character(my_colors['human_HBM'])
rc[grep('human.*CSHL', sample_names)] = as.character(my_colors['human_CSHL'])

par(oma=c(5,5,5,60),xpd=NA)
heatmap.2(common_unmelt, distfun=function(x) dist(x,method='euclidean'),hclustfun=function(x) hclust(x,method='complete'),trace='none',cexRow=2,
	col=rev(heat.colors(16)),keysize=1,labCol='',RowSideColors=rc,lhei=c(1,6), lwid=c(2,6))
legend(x=2.1, y=1,legend=names(my_colors),fill=as.character(my_colors),cex=2)
title(main='Common tissue-specific genes')

sample_names = colnames(ALL_real_thr_expr[,-(1:2)])
rc = rep(0,length(sample_names))
rc[grep('mouse', sample_names)] = as.character(my_colors['mouse'])
rc[grep('HBM', sample_names)] = as.character(my_colors['human_HBM'])
rc[grep('human.*CSHL', sample_names)] = as.character(my_colors['human_CSHL'])

#par(oma=c(5,5,5,60),xpd=NA)
#heatmap.2(as.matrix(log2(t(ALL_real_thr_expr[,-(1:2)]))), distfun=function(x) dist(x,method='euclidean'), 
#	hclustfun=function(x) hclust(x,method='complete'),trace='none',cexRow=2, Rowv=NULL, Colv=NULL,
#	col=rev(heat.colors(16)),keysize=1,labCol='',lhei=c(1,6), lwid=c(2,6))
#legend(x=2.1, y=1,legend=names(my_colors),fill=as.character(my_colors),cex=2)
#title(main='Common tissue-specific genes')

par(oma=c(5,5,5,60),xpd=NA)
heatmap.2(cor(ALL_real_thr_expr[,-(1:2)], u='p',m='s'), distfun=function(x) dist(x,method='euclidean'),
	hclustfun=function(x) hclust(x,method='complete'),trace='none',cexRow=2,
	col=rev(heat.colors(16)),keysize=1,labCol='',RowSideColors=rc,lhei=c(1,6), lwid=c(2,6))
legend(x=2.1, y=1,legend=names(my_colors),fill=as.character(my_colors),cex=2)
title(main='Common tissue-specific genes')


dev.off()




