# This script is useful for:
# normalizing samples by quantile


##------------
## LIBRARIES
##------------ 

suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)
pseudocount = 1e-05
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-m", "--metadata"), help="tsv file with the metadata"),
make_option(c("-b", "--bottom"), type='double', help="fraction to remove from the bottom sorted distribution column-wise"),
make_option(c("-t", "--top"), type='double', help="fraction to remove from the top sorted distribution column-wise"),
make_option(c("-o", "--output"), help="output file name WIHTOUT extension [default=%default]", default="qnorm")
#make_option(c("-s", "--scale_by"), help="choose one or multiple attributes you want to scale by"),
#make_option(c("-r", "--row_first"), action="store_true", help="scale first by rows then by columns", default=FALSE),
#make_option(c("-n", "--n_iter"), type='integer', help="how many times to iterate [default=%default]", default=20)
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)

###############
# BEGIN
###############

# read options
m <- read.table(opt$input_matrix, h=T)
mdata <- read.table(opt$metadata, h=T, row.names=NULL, sep="\t")
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)

if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}


# apply the log if required
#if (opt$log) { m = log2(replace(m, is.na(m), 0) + pseudocount) }

gp1 = ggplot(melt(m), aes(x=value)) + geom_density() + facet_wrap(~variable)

# sort the columns of m, be careful to keep the labels
# I will create two matrices, one with the sorted values, the other with the sorted labels

sort_m = apply(m, 2, sort, na.last=T, d=F)
sort_genes_m = apply(m, 2, function(x) as.numeric(rownames(sort_df(as.data.frame(x)))))


shift1 = opt$bottom
shift2 = opt$top

# from this extract a range 
sort_m_range = sort_m[floor(nrow(sort_m)*(shift1)) : floor(nrow(sort_m)*(shift2)),]
sort_genes_m_range = sort_genes_m[floor(nrow(sort_genes_m)*(shift1)) : floor(nrow(sort_genes_m)*(shift2)),]

gp2 = ggplot(melt(as.data.frame(sort_m_range)), aes(x=value)) + geom_density() + facet_wrap(~variable)

# calculate mean and standard deviation for the gated matrix
#scl_mean = apply(sort_m_range, 2, mean) 
#scl_sd = apply(sort_m_range, 2, sd)

# with mean and standard deviation scale the whole sorted and plot
#scl_sort_m = sapply(1:ncol(sort_m), function(i) (sort_m[,i]-scl_mean[i])/scl_sd[i]) 
#colnames(scl_sort_m) <- colnames(sort_m)

scl_sort_m = sort_m_range

#gp3 = ggplot(melt(as.data.frame(scl_sort_m)), aes(x=value)) + geom_density() + facet_wrap(~variable)

# do the average by rows and reassign the gene labels
av = apply(scl_sort_m, 1, mean)
qnorm_m = apply(sort_genes_m_range, 2, function(col) av[order(col)])
gp4 = ggplot(melt(as.data.frame(qnorm_m)), aes(x=value)) + geom_density() + facet_wrap(~variable)

# write the output matrix to a file
#if (length(char_cols) != 0) {qnorm_m <- cbind(genes, qnorm_m)}
write.table(qnorm_m, sprintf("%s.tsv", opt$output), quote=F, sep='\t', row.names=F)

# plot all the distributions after each of the steps
pdf(sprintf("%s.pdf", opt$output))
gp1;gp2;#gp3;
gp4
dev.off()

q(save='no')
