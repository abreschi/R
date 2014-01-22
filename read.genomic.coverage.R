#!/usr/bin/env Rscript

#options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 



##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze. It must have a header."),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-M", "--merge_mdata_on"), help="which field corresponds to the ids in the summary file [default=%default]", default="labExpId"),
make_option(c("-o", "--output"), help="output file name (without extension) [default=%default]", default="summary.out"),
make_option(c("-f", "--facet"), help="dashboard field by which the individuals are faceted")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)


##------------
## LIBRARIES
##------------ 

cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
cat("DONE\n\n")



##--------------------##
## BEGIN              ##
##--------------------##

# read the matrix from the command line
m = read.table(opt$input_matrix, h=F, col.names=c("labExpId", "type", "region", "nb_reads"))

# read the metadata from the metadata file
if (!is.null(opt$metadata)) {
	mdata = read.table(opt$metadata, h=T, sep='\t')
	mdata[opt$merge_mdata_on] <- sapply(mdata[opt$merge_mdata_on], function(x) gsub(",", ".", x))
}

# separate total from the rest
df = merge(m, setNames(subset(m, type=="total")[c(1,4)], c("labExpId", "total")), by="labExpId")

# merge split and continuous
all = merge(aggregate(nb_reads~labExpId+region, subset(df, region!="total"), sum), subset(m, type=="total")[c(1,2,4)], by="labExpId")
colnames(all)[c(3,5)] <- c("nb_reads", "total")
all$type <- "all"
df = rbind(df, all)

# attach the metadata
if (!is.null(opt$metadata)) {
	mdata_header = unique(c(opt$facet, opt$merge_mdata_on))
	df = merge(df, unique(mdata[mdata_header]), by.x='labExpId', by.y=opt$merge_mdata_on)
}


# ----------------- ggplot options ------------------------------

theme_set(theme_bw(base_size=16))

gp = ggplot(subset(df, type!="total"), aes(y=nb_reads/total*100, x=region))
gp = gp + geom_boxplot(aes(color=type), fill=grey)
if(!is.null(opt$metadata)) {gp = gp + facet_grid(as.formula(sprintf("~%s", opt$facet)))}
gp = gp + labs(y='Proportion of mapped reads (%)', x="")
gp = gp + theme(axis.text = element_text(size=13, angle=45, h=1))
gp = gp + scale_color_brewer(palette="Set1")

w=8; h=6

ggsave(filename=sprintf("%s.pdf", opt$output), h=h, w=w)
ggsave(filename=sprintf("%s.png", opt$output), h=h, w=w)
ggsave(filename=sprintf("%s.eps", opt$output), h=h, w=w)


q(save='no')
