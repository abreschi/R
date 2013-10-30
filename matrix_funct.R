
# This script is useful for:
# normalizing samples by row scaling


##------------
## LIBRARIES
##------------
cat("NOTE: Treat Inf and -Inf as NAs\n\n")
 
cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
source("~/R/functions.R")
cat("DONE\n\n")


options(stringsAsFactors=F)
pseudocount = 1e-05
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


# ==================
# DEBUG OPTIONS
# ==================


opt = list()
opt$input_matrix = "/users/rg/abreschi/Documents/human-mouse/antisense-transcription/s_as_read_counts/corr-mean/log2_FALSE.pdcn_1e-04.NAs_TRUE.human.mouse.S_AS_ratio.tsv"
opt$metadata = '/users/rg/abreschi/Documents/human-mouse/paper-sample-clustering/merged_RNA_dashboard_files.crg.tsv'
opt$mean_by="organism"
opt$output_suffix = "s_as_ratio"
opt$func = "mean"
opt$not_na = 0.7
opt$log=FALSE
opt$replace_na = FALSE

##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-k", "--replace_na"), action="store_true", default=FALSE, help="use this if you want NAs to be replaces by 0 [default=%default]"),
make_option(c("-p", "--pseudocount"), type="double", help="specify a pseudocount for the log [default=%default]", default=pseudocount),
make_option(c("-m", "--metadata"), help="tsv file with the metadata"),
make_option(c("-s", "--mean_by"), help="choose one or multiple attributes you want to average by"),
make_option(c("-o", "--output_suffix"), help="additional output suffix [default=%default]", default='out'),
make_option(c("-f", "--func"), help="choose the function <mean>, <sd>, <sum>, <median>, <entropy> [default=%default]", default='mean'),
make_option(c("-n", "--not_na"), help="fraction of not NA values in the vector for the mean [default=%default]", default=1, type='double')
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
if (opt$replace_na) {m<-replace(m, is.na(m), 0)}
if (opt$log) { m = log10(replace(m, is.na(m), 0) + opt$pseudocount) }

mdata = subset(mdata, labExpId %in% colnames(m))

# ----------------
# Functions
# ----------------

my_mean = function(x) {
ifelse((sum(!is.na(x)) < (opt$not_na*length(x))), NA, mean(x,na.rm=T) )
}

my_sd = function(x) {
ifelse((sum(!is.na(x)) < (opt$not_na*length(x))), NA, sd(x,na.rm=T) )
}

my_sum = function(x) {
ifelse((sum(!is.na(x)) < (opt$not_na*length(x))), NA, sum(x,na.rm=T) )
}

my_median = function(x) {
ifelse((sum(!is.na(x)) < (opt$not_na*length(x))), NA, median(x,na.rm=T) )
}

my_entropy = function(x) {
ifelse((sum(!is.na(x)) < (opt$not_na*length(x))), NA, entropy(na.omit(x)) )
}


if (opt$func == "mean") {func = my_mean}
if (opt$func == "sd") {func = my_sd}
if (opt$func == "sum") {func = my_sum}
if (opt$func == "median") {func = my_median}
if (opt$func == "entropy") {func = my_entropy}


# apply the function to the whole matrix if no value is provided
if (is.null(opt$mean_by)) {
new_m = apply(m, 1, func); colnames(new_m) <- c(opt$func)
} else {
# apply the function to the levels of the specified factors
mean_by = strsplit(opt$mean_by, ",")[[1]]
meltm = melt(as.matrix(m), varnames = c("gene_index", "labExpId"))
df = merge(meltm, unique(mdata[c("labExpId",mean_by)]), by = "labExpId")
df$value[abs(df$value)==Inf] <- NA
aggr = aggregate(as.formula(sprintf("value~gene_index+%s", paste(mean_by,collapse="+"))), df, func, na.action="na.pass")
aggr = dcast(aggr, as.formula(sprintf("gene_index~%s", paste(mean_by,collapse="+"))))
if (length(char_cols)==0) {new_m = aggr; colnames(new_m)[1] = "gene"}
if (length(char_cols)!=0) {new_m = merge(genes, aggr, by.y="gene_index", by.x="row.names")[,-1]}
}



#--------------
# print output
#--------------

output = sprintf('%s.by_%s.n_%s.%s', opt$func, opt$mean_by, opt$not_na, opt$output_suffix)
write.table(new_m, sprintf('%s.tsv',output), quote=F, sep='\t', row.names=F)
#pdf(sprintf("%s.pdf", output)); gp1; gp2; dev.off()
q(save='no')
