#!/usr/bin/env Rscript 

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin", 
	help="the matrix with READ COUNTS you want to analyze. Should have a header (n-1 columns) [default=%default]"),

make_option(c("-M", "--method"), type="character", help="Normalization method [default=%default]
	
		cpm      :
		rlog     :
		divsum   : divide by the column totals (useful to convert from RPKM to TPM)
		len      : divide by the fragment length * 10^3 (pipe the output to cpm to get RPKM)
		quantile : quantile normalization (normalize.quantiles() function from package preprocessCore)
"),

make_option(c("-s", "--scaling_factors"), default="TMM", help="How to compute scaling factors [default=%default]
	
		TMM    :
		none   :
"),

make_option(c("-p", "--pseudocount"), type="double", help="A pseudocount to add before normalizing with cpm"),

make_option(c("-L", "--lengths"), help="Two-column file with no header. 
		col1: gene id (same as matrix rows), col2: length [default=%default]"),
#make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
#make_option(c("-G", "--merge_mdata_on"), default="labExpId",
#	help="Column in the metadata with the header of the input matrix [default=%default]"),
#make_option(c("-f", "--formula"), help="formula"),
make_option(c("-t", "--total"), type="integer", help="Filter by total count per gene > t [default=%default]"),

make_option(c("-S", "--lib_sizes"), help="Two-column file with no header. col1: header of matrix, col2: library sizes.
		 If not provided, the sum of the column will be used as library size"),

make_option(c("-N", "--output.norm"), help="File name for normalization factors"),
make_option(c("-o", "--output"), default="stdout", help="output file name [default=%default]"),
make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="verbose output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list, description="\nNormalize a matrix")
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

suppressPackageStartupMessages(library("methods"))
suppressPackageStartupMessages(library("data.table"))


##--------------------##
## BEGIN              ##
##--------------------##

merge_mdata_on = opt$merge_mdata_on

# read the matrix from the command line
#if(opt$input_matrix == "stdin"){inF='file:///dev/stdin'}else{inF=opt$input_matrix}
#m = data.frame(fread(inF, header=T), row.names=1)
if(opt$input_matrix == "stdin"){inF=file("stdin")}else{inF=opt$input_matrix}
m = read.table(inF, h=T, sep="\t")

# Replace missing values with 0
m = replace(m, is.na(m), 0)

# Filter by total number of reads per gene if asked
if (!is.null(opt$total)) {
	m = m[rowSums(m)>opt$total, ]
}

## =========================== Metadata =======================
#
#if (!is.null(opt$metadata)) {
#	# read the metadata
#	mdata = read.table(opt$metadata, h=T, sep="\t")
#	# Get the fields from the formula 
#	fields = strsplit(sub("~", "", opt$formula), split="[+:*]")[[1]]
#	merge_mdata_on = "labExpId"
#	# Format the metadata
#	mdata = unique(mdata[unique(c(merge_mdata_on, fields))])
#	rownames(mdata) <- mdata[,merge_mdata_on]
#	mdata <- mdata[match(colnames(m), mdata[,merge_mdata_on]), drop=FALSE]
#	design = model.matrix(eval(as.formula(opt$formula)), data=mdata)
#}


# ========================= Gene lengths =======================

if (opt$method == "len") {
	if (is.null(opt$lengths)) {
		cat("ERROR: Please provide gene lengths\n")
		q(save='no')
	}
	geneLen = read.table(opt$lengths, h=F, sep="\t")
	if (!all(rownames(m) %in% geneLen[,1])) {
		cat("ERROR: Some elements do not have a length\n")
		q(save='no')
	}
	geneLen = geneLen[geneLen[,1] %in% rownames(m),]
	geneLen = geneLen[match(rownames(m), geneLen[,1]),]
	out = (m/geneLen[,2])*1000
}


# ======================= Scaling factors =======================

if (opt$method %in% c("cpm")) {

	# Load edgeR
	suppressPackageStartupMessages(library(edgeR))

	# Convert all the values of the matrix to integer (because we want counts)
	m[1:ncol(m)] <- apply(m, 2, as.integer)
	# Create count object for edgeR
	M = DGEList(m)

	# Check for user-provided library sizes
	if (!is.null(opt$lib_sizes)) {
		lib.sizes = read.table(opt$lib_sizes, h=F, sep="\t")
		lib.sizes = lib.sizes[match(lib.sizes$V1, colnames(m)), "V2"]
		M$samples$lib.size <- lib.sizes
	}

	# ****************
	#      TMM       
	# ****************     
	
	if (opt$scaling == "TMM") {
		M <- calcNormFactors(M, method="TMM")
		if (!is.null(opt$output.norm)) {
			normFactors = data.frame(a=colnames(m), b=M$samples$norm.factors)
			write.table(normFactors, file=opt$output.norm, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
		}
	}
	
	# ****************
	#      none       
	# ****************     

	if (opt$scaling == "none") {
		M$samples$norm.factors <- rep(1, ncol(m))
	}

}

# ======================= Normalization =======================

# ****************
#      divsum       
# ****************     

if (opt$method == "divsum") {
	out = sweep(m, 2, colSums(m), FUN="/")*1e+06
}


# ****************
#      cpm       
# ****************     

if (opt$method == "cpm") {
	if (!is.null(opt$pseudocount)) {
		M$counts <- M$counts + opt$pseudocount
	}
#	out <- cpm(M, normalized.lib.sizes=TRUE)
	out <- sweep(M$counts, 2, M$samples$lib.size * M$samples$norm.factors, FUN="/") * 1e+06
}


# # *********************
# #      voom 
# # *********************
# 
# if (opt$method == "voom") {
# 
# 	# Load limma
# 	suppressPackageStartupMessages(library(limma))
# 
# 	if (is.null(opt$metadata)) {
# 		mdata = data.frame(a = colnames(m), row.names = colnames(m))
# 		design <- model.matrix(~1, data=mdata)
# 	} 
# 	
# 	v <- voom(M, design, plot=FALSE)
# 	out <- v$E
# }
# 

# ******************
#     rlog
# ******************

if (opt$method == "rlog") {

	# Convert all the values of the matrix to integer (because we want counts)
	m[1:ncol(m)] <- apply(m, 2, as.integer)

#	# read the metadata from the metadata file
#	mdata = read.table(opt$metadata, h=T, sep='\t')
#
#	# specify the design to the program
#	fields = strsplit(opt$fields, ",")[[1]]
#
#
#	if (length(fields) == 1) {
#		mdata = unique(mdata[, c(merge_mdata_on, fields)])
#		colData = mdata
#		rownames(colData) <- colData[, merge_mdata_on]
#		colData = colData[match(colnames(m), rownames(colData)),]
#	condition = factor(sapply(colnames(m), function(x) unique(subset(mdata, labExpId == x)[,opt$fields])))
#	}else{print('cannot handle multiple fields yet');q(save='no')}
	
	# Load DESeq2
	if (opt$verbose) {cat("Loading library... ")}
	suppressPackageStartupMessages(library('DESeq2'))
	if (opt$verbose) {cat("DONE\n")}
		
	# create count object for DESeq
	colData = data.frame(colnames(m), row.names=colnames(m))
	dds = DESeqDataSetFromMatrix(countData = m, colData = colData, design=~1)

	# rlog
	rld <- rlog(dds)
	out <- assay(rld)
}

# *******************
#    quantile
# *******************

if (opt$method == "quantile") {
	suppressPackageStartupMessages(library('preprocessCore'))
	out = normalize.quantiles(as.matrix(m))
	colnames(out) <- colnames(m)
	rownames(out) <- rownames(m)
}


# =================== OUTPUT ======================

out = round(out, digits=5)
outF = ifelse(opt$output=="stdout", "", opt$output)
write.table(out, file=outF, quote=FALSE, sep="\t")

q(save='no')
