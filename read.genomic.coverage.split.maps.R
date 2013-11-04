#!/usr/bin/env Rscript

#options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 



##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
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
m = read.table(opt$input_matrix, h=F, col.names=c("labExpId","n_det_el","element"))

# read the metadata from the metadata file
if (!is.null(opt$metadata)) {
	mdata = read.table(opt$metadata, h=T, sep='\t')
	mdata[opt$merge_mdata_on] <- sapply(mdata[opt$merge_mdata_on], function(x) gsub(",", ".", x))
}

# sum split maps and exonic
df_ex = merge(subset(m, element=="exonic_reads"), subset(m, element=="split_reads"), by="labExpId")
df_ex$n_det_el.x = with(df_ex, n_det_el.x+n_det_el.y)
m[which(m$element=="exonic_reads"),] <- df_ex[1:3]
m = m[-which(m$element=="split_reads"),]

# separate total from the rest
mtot = subset(m, element=='total_reads')
m = subset(m , element!='total_reads')

# attach the total as additional column
df = merge(m, mtot, by='labExpId')

# attach the metadata
mdata_header = unique(c(opt$facet, opt$merge_mdata_on))
df = merge(df, unique(mdata[mdata_header]), by.x='labExpId', by.y=opt$merge_mdata_on)
#df$element.x <- factor(df$element.x, levels=c('exonic_reads', 'intronic_reads', 'exonic_intronic_reads', 'genic_reads', 'intergenic_reads'))
df$element.x <- gsub("_reads", "", df$element.x)
df$element.x <- factor(df$element.x, levels=c('exonic', 'intronic', 'exonic_intronic', 'genic', 'intergenic'))


# ----------------- ggplot options ------------------------------

theme_set(theme_bw(base_size=16))

gp = ggplot(df) 
gp = gp + geom_boxplot(aes(y=n_det_el.x/n_det_el.y*100, x=element.x), fill="grey") 
gp = gp + facet_grid(as.formula(sprintf("~%s", opt$facet)))
gp = gp + labs(y='Proportion of mapped reads (%)', x="")
gp = gp + theme(axis.text = element_text(size=13, angle=45, h=1))

w=5
h=5

ggsave(filename=sprintf("%s.pdf", opt$output), h=h, w=w)
ggsave(filename=sprintf("%s.png", opt$output), h=h, w=w)
ggsave(filename=sprintf("%s.eps", opt$output), h=h, w=w)


q(save='no')
