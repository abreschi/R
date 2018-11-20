#!/usr/bin/env Rscript

cat("NOTE: Treat Inf and -Inf as NAs\n\n", file=stderr())
 
suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)


# ==================
# DEBUG OPTIONS
# ==================

opt = list()
opt$input_matrix = "/users/rg/abreschi/Documents/human-mouse/antisense-transcription/s_as_read_counts/corr-mean/log2_FALSE.pdcn_1e-04.NAs_TRUE.human.mouse.S_AS_ratio.tsv"
opt$metadata = '/users/rg/abreschi/Documents/human-mouse/paper-sample-clustering/merged_RNA_dashboard_files.crg.tsv'
opt$mean_by="organism"
opt$output = "s_as_ratio"
opt$func = "mean"
opt$not_na = 0.7
opt$log=FALSE
opt$replace_na = FALSE

##################
# OPTION PARSING
##################

option_list <- list(

make_option(c("-i", "--input_matrix"), default="stdin",
	help="the matrix you want to analyze. Stdin to read from stdin"),

make_option(c("--dt"), action="store_true", default=F,
	help="read the matrix as data.table"),

make_option(c("-l", "--log"), action="store_true", default=FALSE, 
	help="apply the log10, before applying the function [default=%default]"),

make_option(c("-k", "--replace_na"), action="store_true", default=FALSE, 
	help="use this if you want NAs to be replaced by 0 [default=%default]"),

make_option(c("-p", "--pseudocount"), type="double", default=1e-04,
	help="specify a pseudocount for the log [default=%default]"),

make_option(c("-m", "--metadata"), 
	help="tsv file with the metadata"),

make_option(c("-s", "--mean_by"), 
	help="choose one or multiple attributes you want to average by"),

make_option(c("-o", "--output"), default="out.tsv", 
	help="Output file name. stdout for printing on standard output [default=%default]"),

make_option(c("-f", "--func"), default="mean",
	help="choose the function <mean>, <sd>, <sum>, <median>, <entropy>, <nentropy>, <tau>, <cv> [default=%default]"),

make_option(c("-C", "--byColumns"), action="store_true", default=FALSE,
	help="apply the function to the columns, instead of rows [default=%default]"),

make_option(c("-n", "--not_na"), type="double", default=1, 
	help="fraction of not NA values in the vector for the mean. If NAs are replaced they are not counted [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="verbose output")
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list,
	description = "\nExecute a function on matrix rows or columns"
	)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


##------------
## LIBRARIES
##------------

if (opt$verbose) {cat("Loading libraries... ", file=stderr())}
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
source("./functions.R")
if (opt$verbose) {cat("DONE\n\n", file=stderr())}


# ==========================================
# Function for loading Rdata
# ==========================================

load_obj <- function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}

###############
# BEGIN
###############

# read table
if (opt$verbose) {cat(sprintf("%s: ", Sys.time()), "Reading matrix... ")}
if (opt$dt) {
	inF = ifelse(opt$input_matrix == "stdin", "file:///dev/stdin", opt$input_matrix)
	m = fread(inF)
	row.names = m[,1][[1]]
	m = m[,-1]
} else {
	if (opt$input_matrix == "stdin") {
	    m = read.table(file("stdin"), h=T)
	} else {
		m <- try(load_obj(opt$input_matrix), silent=T)
	    if (class(m) == "try-error") {m <- read.table(opt$input_matrix)}
	}
	row.names = rownames(m)
}
if (opt$verbose) {cat("DONE\n")}

# apply the log if required
if (opt$replace_na) {m <- replace(m, is.na(m), 0)}
if (opt$log) {m = log10(m + opt$pseudocount)}

# Set the result to NA when too many missing values are present
func = function(x) {
	ifelse((sum(!is.na(x)) < (opt$not_na*length(x))), NA, format(eval(parse(text=opt$func))(x,na.rm=T), digits=5))
}


# Apply the function by columns
if (opt$byColumns) {
	if (!is.null(opt$metadata)) {
		mdata <- read.table(opt$metadata, h=T, sep="\t", quote="", check.names=FALSE)
	}
	if (!is.null(opt$mean_by)) {
		mean_by = strsplit(opt$mean_by, ",")[[1]]
		m[,ncol(m)+1] <- mdata[,mean_by][match(rownames(m), mdata[,"gene"])]
		colnames(m)[ncol(m)] <- mean_by
		form = as.formula(sprintf(".~%s", mean_by))
		new_m <- aggregate(form, m, func)
	} else {
		new_m = data.frame(id=colnames(m), fun=apply(m, 2, func))
		
	}

} else {

	# apply the function to the whole matrix if no value is provided
	if (is.null(opt$mean_by)) {
		new_m = setNames(data.table("id"=row.names, "func"=apply(m, 1, func)), c("id", opt$func))
	} else {
	
		# apply the function to the levels of the specified factors
		mean_by = strsplit(opt$mean_by, ",")[[1]]
		df = melt(as.matrix(m), varnames = c("gene_index", "labExpId"))
		
		# read metadata and merge with data.frame if needed
		if (!is.null(opt$metadata)) {
			mdata <- read.table(opt$metadata, h=T, row.names=NULL, sep="\t", quote="", check.names=FALSE)
			mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub("[-,]", ".", x))
			mdata = subset(mdata, labExpId %in% colnames(m))
			df = merge(df, unique(mdata[c("labExpId", mean_by)]), by = "labExpId")
		}
		
		df$value[abs(df$value)==Inf] <- NA
		aggr = aggregate(as.formula(sprintf("value~gene_index+`%s`", paste(mean_by,collapse="+"))), df, func, na.action="na.pass")
		aggr = dcast(aggr, as.formula(sprintf("gene_index~`%s`", paste(mean_by,collapse="+"))))
		new_m = aggr
	
	#	if (length(char_cols)==0) {new_m = cbind(gene=genes, new_m)}
	#	if (length(char_cols)!=0) {new_m = merge(genes, new_m, by.y="gene_index", by.x="row.names")[,-1]}
	}
}	


#--------------
# print output
#--------------

output = ifelse(opt$output == "stdout", "", opt$output)
write.table(new_m, file = output, quote=F, sep='\t', row.names=F)

q(save='no')

