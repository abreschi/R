
# This script is useful for:
# normalizing samples by row scaling


##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq))
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
make_option(c("-s", "--scale_by"), help="choose one or multiple attributes you want to scale by. Null for normalize the whole matrix")
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

m = (apply(m, 2, as.integer))

# apply the log if required
if (opt$log) { m = log2(replace(m, is.na(m), 0) + pseudocount) }

mdata = subset(mdata, labExpId %in% colnames(m))



# store the plot of the distribution with input matrix
gp1 = ggplot(melt(m), aes(x=value)) + geom_density() + facet_wrap(~variable)



##########
# EQUILIBRATION
##########


DESeq_norm = function(matr) {
	condition = rep('group',ncol(matr))
	cds = newCountDataSet(matr, condition)
	cds = estimateSizeFactors(cds)
	return(counts(cds, normalized=T))
}

new_m = m

# scale the whole matrix if no value is provided
if (is.null(opt$scale_by)) {
new_m = DESeq_norm(new_m)
}

# scale the sub-matrices defined the scale_by option
if (!is.null(opt$scale_by)) {
scale_by <- strsplit(opt$scale_by, ",")[[1]]
if (length(scale_by) != 1){
ids = apply(data.frame(unique(mdata[, scale_by])), 1, function(x) unique(merge(t(as.data.frame(x)), mdata, by=scale_by)$labExpId ))}
if (length(scale_by) == 1){
ids = sapply(unique(mdata[, scale_by]), function(x) unique(mdata[ mdata[,scale_by] == x,]$labExpId))}
# apply normalization
if (length(scale_by) != 1){for (i in 1:length(ids)) { new_m[, ids[[i]]] <- DESeq_norm(new_m[,ids[[i]]])}  }
if (length(scale_by) == 1){for (i in 1:ncol(ids)) { new_m[, ids[[i]]] <- DESeq_norm(new_m[,ids[[i]]])}  }
}


if (is.null(opt$scale_by)) {gp2 = ggplot(melt(new_m), aes(x=value)) + geom_density() + facet_wrap(~Var2)}else{
gp2 = ggplot(melt(new_m), aes(x=value)) + geom_density() + facet_wrap(~variable)}
if (length(char_cols) != 0) {new_m <- cbind(genes, new_m)}


# print output
#--------------
output = sprintf('DESeq_norm.log_%s.pscnt_%s', opt$log, opt$pseudocount)
write.table(new_m, sprintf('%s.tsv',output), quote=F, sep='\t', row.names=F)
pdf(sprintf("%s.pdf", output)); gp1; gp2; dev.off()
q(save='no')
