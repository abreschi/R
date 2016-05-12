#!/usr/bin/env Rscript

##------------
## LIBRARIES
##------------ 
cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(plyr))
cat("DONE\n\n")

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-A", "--input1"), help="the first matrix of samples"),
make_option(c("-B", "--input2"), help="the second matrix of samples, columns must have same order as in input1"),
make_option(c("--tag1"), help="short tag describing the first input"),
make_option(c("--tag2"), help="short tag describing the second input"),
make_option(c("-L", "--labels"), help="select the labels from the dashboard fields, comma-separated"),
make_option(c("-l", "--log"), help="apply the log to <both>, <A>, <B>, <none> [default=%default]", default="none"),
make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-m", "--metadata"), default=NULL, help="tsv file with metadata on matrix experiment"),
make_option(c("-r", "--representation"), help="choose <scatter>, <hex> [default=%default]", default="hex"),
make_option(c("-n", "--replace_na_with_0"), action="store_true", help="replace NAs with 0 [default=%default]", default=FALSE),
make_option(c("--facet_nrow"), default=NULL, type="integer", help="Number of rows for faceting [default=%default]"),
make_option(c("-o", "--output"), help="output file name without extension", default='out'),
make_option(c("-v", "--verbose"), default=FALSE, action="store_true", help="verbose")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {
	print(opt)
}


##--------------------##
## BEGIN              ##
##--------------------##

#output = sprintf("log_%s.pseudo_%s.fillby_%s.%s", opt$log, opt$pseudocount, opt$fill_by, opt$output)

m1 = read.table(opt$input1, h=T)
m2 = read.table(opt$input2, h=T)

if (opt$replace_na_with_0) {
m1 = replace(m1, is.na(m1), 0)
m2 = replace(m2, is.na(m2), 0)
}

# check order of headers
if (!all(colnames(m1) == colnames(m2))) {print('Headers are not the same!'); q(save='no') }


mm1 = melt(cbind(m1, gene=rownames(m1)), value.name="value1", variable.name='labExpId')
mm2 = melt(cbind(m2, gene=rownames(m2)), value.name="value2", variable.name='labExpId')

df = merge(mm1, mm2, by=c('gene', 'labExpId'))

if (opt$verbose) {print(head(df))}


# Log-transformation
if (opt$log=="both") {df$value1 <- log10(df$value1+opt$pseudocount); df$value2 <- log10(df$value2+opt$pseudocount)}
if (opt$log=="A") {df$value1 <- log10(df$value1+opt$pseudocount)}
if (opt$log=="B") {df$value2 <- log10(df$value2+opt$pseudocount)}


# read the metadata
if (!is.null(opt$metadata)) {
	mdata = read.table(opt$metadata, h=T, sep='\t', row.names=NULL)
	# select only the requested ones
	labels = strsplit(opt$labels, ',')[[1]]
	mdata = unique(mdata[c('labExpId', labels)])
	mdata$labels = apply(mdata[labels], 1, function(x) paste(x, collapse='\n'))
	mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub("[,-]", ".", x))
	# add this info to the data.frame
	df = merge(df, mdata[c('labExpId',"labels")], by='labExpId')

# calculate correlation coefficient
df = ddply(df, .(labExpId,labels), transform, pearson = cor(value1, value2, m='p', use='p'))
} else {
# calculate correlation coefficient
df = ddply(df, .(labExpId), transform, pearson = cor(value1, value2, m='p', use='p'))
df$labels = df$labExpId
}


##########
# PLOTTING
##########

gp = ggplot(df, aes(x=value1, y=value2))
if (opt$representation == "hex") {
gp = gp + geom_hex(aes(fill=cut(..count..,c(0,1,2,5,10,25,50,75,100,500,Inf))),binwidth=c(.2,.2))
gp = gp + scale_fill_manual('counts',values=terrain.colors(10))}
if (opt$representation == "scatter") {
gp = gp + geom_point(shape=".")}
gp = gp + facet_wrap(~labels, nrow=opt$facet_nrow)
gp = gp + geom_abline(slope=1, intercept=0, col='blue')
gp = gp + geom_text(aes(x=min(value1,na.rm=T), y=max(value2,na.rm=T)-0.5, label=sprintf("r=%0.3f",pearson)), color='blue', hjust=0)
gp = gp + labs(x=opt$tag1, y=opt$tag2)
#gp
ggsave(filename=sprintf('%s.pdf', opt$output), h=7, w=7)
ggsave(filename=sprintf('%s.jpg', opt$output), h=7, w=7)



q(save='no')



