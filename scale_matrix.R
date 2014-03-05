#!/usr/bin/env Rscript


# This script is useful for:
# normalizing samples by row scaling

##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin", 
	help="the matrix you want to analyze. \"stdin\" for reading from standard input [default=%default]"),

make_option(c("-l", "--log10"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-p", "--pseudocount"), type="double", help="specify a pseudocount for the log [default=%default]", default=1e-04),
make_option(c("-m", "--metadata"), help="tsv file with the metadata"),
make_option(c("-s", "--scale_by"), help="choose one or multiple attributes you want to scale by"),
make_option(c("-r", "--row_first"), action="store_true", help="scale first by rows then by columns", default=FALSE),
make_option(c("-C", "--center"), action="store_true", help="subtract the mean [default=%default]", default=FALSE),
make_option(c("-S", "--scale"), action="store_true", help="divide by the standard deviation [default=%default]", default=FALSE),
make_option(c("-n", "--n_iter"), type='integer', help="how many times to iterate. Choose 0 for one-dimension scaling [default=%default]", default=20),

make_option(c("-o", "--output"), default="stdout", 
	help="output file name. \"stdout\" to redirect to stdout. [default=%default]"),

make_option(c("--verbose"), action="store_true", default=FALSE,
	help="verbose output. [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


###############
# BEGIN
###############

# read options
if (opt$input_matrix == "stdin") {
	m <- read.table(file("stdin"), h=T)
} else {
	m <- read.table(opt$input_matrix, h=T)
}

# Remove character columns
char_cols <- which(sapply(m, class) == 'character')
if (opt$verbose) {sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)}
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# Check the scale attribute
#if (opt$center_only) {scale_attr = FALSE}

# apply the log if required
if (opt$log10) { m = log10(replace(m, is.na(m), 0) + opt$pseudocount) }


# Read the metadata
if (!is.null(opt$metadata)) {
mdata <- read.table(opt$metadata, h=T, row.names=NULL, sep="\t")
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))
mdata = subset(mdata, labExpId %in% colnames(m))
}




# store the plot of the distribution with input matrix
#gp1 = ggplot(melt(m), aes(x=value)) + geom_density() + facet_wrap(~variable)



##########
# EQUILIBRATION
##########

# Comments on scale() function from R:
# - if the vector has sd==0, all the scaled values will be NAs.
# - if the vector contains NAs, scale omits them.

# scaling by rows and by columns in this order for 20 times

if (opt$row_first) {
# scale first by rows then by columns iteratively
equil = function(matr) {
	if (opt$n_iter == 0) {matr = t(scale(t(matr), center=opt$center, scale=opt$scale))
		}else{
		for (i in 1:opt$n_iter) {
			matr = t(scale(t(matr), center=opt$center, scale=opt$scale))
			matr = scale(matr, center=opt$center, scale=opt$scale)
			}
		}
	return(matr)
	}
}else{
# scale first by columns then by rows iteratively
equil = function(matr) {
	if (opt$n_iter == 0) {matr = scale(matr, center=opt$center, scale=opt$scale)
		}else{
		for (i in 1:opt$n_iter) {
			matr = scale(matr, center=opt$center, scale=opt$scale)
			matr = t(scale(t(matr), center=opt$center, scale=opt$scale))
			}
		}
	return(matr)
	}
}

new_m = m

# scale the whole matrix if no value is provided
if (is.null(opt$scale_by)) {
new_m = equil(new_m)
}

# scale the sub-matrices defined the scale_by option
if (!is.null(opt$scale_by)) {
scale_by <- strsplit(opt$scale_by, ",")[[1]]
if (length(scale_by) != 1){
ids = apply(data.frame(unique(mdata[, scale_by])), 1, function(x) unique(merge(t(as.data.frame(x)), mdata, by=scale_by)$labExpId ))}
if (length(scale_by) == 1){
ids = sapply(unique(mdata[, scale_by]), function(x) unique(mdata[ mdata[,scale_by] == x,]$labExpId))}
# apply normalization
if (length(scale_by) != 1){for (i in 1:length(ids)) { new_m[, ids[[i]]] <- equil(new_m[,ids[[i]]])}  }
if (length(scale_by) == 1){for (i in 1:ncol(ids)) { new_m[, ids[[i]]] <- equil(new_m[,ids[[i]]])}  }
}


#if (is.null(opt$scale_by)) {gp2 = ggplot(melt(new_m), aes(x=value)) + geom_density() + facet_wrap(~Var2)}else{
#gp2 = ggplot(melt(new_m), aes(x=value)) + geom_density() + facet_wrap(~variable)}

print_rownames = TRUE
if (length(char_cols) != 0) {new_m <- cbind(genes, new_m); print_rownames=FALSE}


# print output
#--------------
#output = sprintf('equil.log_%s.pscnt_%s.rowFirst_%s', opt$log10, opt$pseudocount, opt$row_first)
#write.table(new_m, sprintf('%s.tsv',output), quote=F, sep='\t', row.names=F)
#pdf(sprintf("%s.pdf", output)); gp1; gp2; dev.off()

if (opt$output == "stdout") {
	output = ""
} else {
	output = opt$output
}


write.table(new_m, output, quote=F, sep="\t", row.names=print_rownames)

q(save='no')

