
# This script is useful for:
# normalizing samples by row scaling


##------------
## LIBRARIES
##------------ 
cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
cat("DONE\n\n")


options(stringsAsFactors=F)
pseudocount = 1e-05
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
#make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-m", "--metadata"), help="tsv file with the metadata"),
make_option(c("-b", "--by"), help="choose ONE factor you want to do the test by"),
make_option(c("-o", "--output_suffix"), help="additional output suffix [default=%default]", default='out'),
#make_option(c("-f", "--func"), help="choose the function <mean>, <sd> [default=%default]", default='mean'),
make_option(c("-n", "--not_na"), help="fraction of not NA values in the vector for the mean [default=%default]", default=1, type='double')
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

mdata = subset(mdata, labExpId %in% colnames(m))

# ----------------
# Functions
# ----------------

#my_mean = function(x) {
#ifelse((sum(!is.na(x)) < (opt$not_na*length(x))), NA, mean(x,na.rm=T) )
#}
#
#my_sd = function(x) {
#ifelse((sum(!is.na(x)) < (opt$not_na*length(x))), NA, sd(x,na.rm=T) )
#}
#
#if (opt$func == "mean") {func = my_mean}
#if (opt$func == "sd") {func = my_sd}


# apply the function to the whole matrix if no value is provided
#if (is.null(opt$mean_by)) {
#new_m = rowMeans(m, na.rm=T); colnames(new_m) <- c("mean")
#}

# scale the sub-matrices defined the scale_by option
if (!is.null(opt$by)) {
by <- strsplit(opt$by, ",")[[1]]
if (length(by) != 1){
#ids = apply(data.frame(unique(mdata[, mean_by])), 1, function(x) unique(merge(t(as.data.frame(x)), mdata, by=mean_by)$labExpId ))
cat("ERROR: Only ONE factor must be provided!\n"); q(save='no')}
if (length(by) == 1){
ids = sapply(unique(mdata[, by]), function(x) unique(mdata[ mdata[,by] == x,]$labExpId))}
# apply normalization
#if (length(mean_by) != 1){cat('to be implemented\n\n');for (i in 1:length(ids)) { new_m[, ids[[i]]] <- equil(new_m[,ids[[i]]])}  }
#if (length(by) == 1){new_m = matrix(ncol=ncol(ids), nrow=nrow(m)); for (i in 1:ncol(ids)) {  
#colnames(new_m) <- colnames(ids); new_m[,i] <- apply(m[,ids[,i]], 1, mean, na.rm=T)}  }
}

wilc_pv = apply(m, 1, function(x) tryCatch({wt = wilcox.test(x[ids[,1]], x[ids[,2]]);wt$p.value}, error = function(e) NA_real_))
ttes_pv = apply(m, 1, function(x) tryCatch({wt = t.test(x[ids[,1]], x[ids[,2]]);wt$p.value}, error = function(e) NA_real_))


#if (is.null(opt$scale_by)) {gp2 = ggplot(melt(new_m), aes(x=value)) + geom_density() + facet_wrap(~Var2)}else{
#gp2 = ggplot(melt(new_m), aes(x=value)) + geom_density() + facet_wrap(~variable)}

#if (length(char_cols) != 0) {new_m <- cbind(genes, new_m)}
df = data.frame(cbind(genes, wilc_pv))


#--------------
# print output
#--------------

output = sprintf('%s.by_%s.n_%s.%s', opt$func, opt$by, opt$not_na, opt$output_suffix)
write.table(new_m, sprintf('%s.tsv',output), quote=F, sep='\t', row.names=F)
#pdf(sprintf("%s.pdf", output)); gp1; gp2; dev.off()
q(save='no')
