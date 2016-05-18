#!/usr/bin/env Rscript 

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin", 
	help="the matrix you want to analyze [default=%default]"),


make_option(c("-l", "--log"), default=NULL, 
	help="The base of the logarithm to tranform the data before batch removal. If NULL no log-transformation is applied [default=%default]"),

make_option(c("-p", "--pseudocount"), type="double", default=0, 
	help="A pseudocount to add when log-transforming [default=%default]"),

make_option(c("-M", "--method"), default="voom", help="The method you want to use for removing the batch effect [default=%default]

		voom   : apply voom to a matrix of read counts to estimate the weigths
			     and then uses the function removeBatchEffect() from limma
		limma  : uses the function removeBatchEffect() from limma on a matrix of 
		         already normalized values
		combat : uses the function combat() from sva on a matrix of already normalized values
"),

make_option(c("-s", "--scaling_factors"), default="TMM", help="How to compute scaling factors, if the method is voom [default=%default]
	
		TMM    :
		none   :
"),
	

make_option(c("--d1"), default="~1", help="Design for voom [default=%default]"),
make_option(c("--d2"), help="Design for removing the batch effect (not including the batch effect)"),
make_option(c("-b", "--batch"), help="Column with the batch info"),

make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-G", "--merge_mdata_on"), default="labExpId",
	help="Column in the metadata with the header of the input matrix [default=%default]"),

make_option(c("-t", "--total"), type="integer", help="Filter by total count per gene > t [default=%default]"),
#make_option(c("-F", "--fields"), help="choose the fields you want to use in the differential expression, comma-separated"),
make_option(c("-S", "--lib.sizes"), help="Two-column file with no header. col1: header of matrix, col2: library sizes"),
make_option(c("-N", "--output.norm"), help="File name for normalization factors"),
make_option(c("-R", "--read_counts"), default=FALSE, action="store_true", help="Output reads counts instead of log2(cpm) [default=%default]"),
make_option(c("-o", "--output"), default="stdout", help="output file name [default=%default]"),
make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="verbose output [default=%default]")
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list, 
	description="\nRemove batch effect from a matrix of read counts, or normalized values"
)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}



# LIBRARIES

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))


##--------------------##
## BEGIN              ##
##--------------------##


# read the matrix from the command line
if(opt$input_matrix == "stdin"){inF=file("stdin")}else{inF=opt$input_matrix}
m = read.table(inF, h=T, sep="\t")

# Replace missing values with 0
m = replace(m, is.na(m), 0)

# Log-transform the values if needed
if (!is.null(opt$log)) {
	base = ifelse(opt$log == "e", exp(1), as.double(opt$log))
	m = log(m+opt$pseudocount, base)
}

inputToBatchRm <- m

# ?TODO: Error if there is attempt to log-transform integer

# =========================== Metadata =======================

merge_mdata_on = opt$merge_mdata_on

# read the metadata
mdata = read.table(opt$metadata, h=T, sep="\t", quote=NULL)

# Get the fields from the formula 
if (is.null(opt$batch)) {
	cat("ERROR: please specify the batch variable\n")
	q(save='no')
}
fields = opt$batch
if (opt$d1 != "~1") {
	fields1 = strsplit(sub("~", "", opt$d1), split="[+:*]")[[1]]
	fields = c(fields, fields1)
}
if (opt$d2 != "~1") {
	fields2 = strsplit(sub("~", "", opt$d2), split="[+:*]")[[1]]
	fields = c(fields, fields2)
}
mdata[opt$merge_mdata_on] <- gsub(",", ".", mdata[,opt$merge_mdata_on])


# Check if all the columns are in the metadata
if (sum(!(colnames(m) %in% mdata[,merge_mdata_on])) >0 ) {
	cat("ERROR: Not all column names in the metadata\n")
	q(save="no")
}

# Format the metadata
mdata = unique(mdata[unique(c(merge_mdata_on, fields))])
rownames(mdata) <- mdata[,merge_mdata_on]
mdata <- mdata[match(colnames(m), mdata[,merge_mdata_on]),, drop=FALSE]
if (opt$verbose) {
	print(mdata)
	print(dim(mdata))
}


# ****************
#    voom+limma       
# ****************     

if (opt$method == "voom") {

	# Filter by total number of reads per gene if asked
	if (!is.null(opt$total)) {
		m = m[rowSums(m)>opt$total, ]
	}

	# Convert all the values of the matrix to integer (because we want counts)
	m[1:ncol(m)] <- apply(m, 2, as.integer)
	# Create count object for edgeR
	M = DGEList(m)
	
	
	# Check for user-provided library sizes
	if (!is.null(opt$lib.sizes)) {
		lib.sizes = read.table(opt$lib.sizes, h=F, sep="\t")
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
	
	
	# ****************
	#      voom
	# ****************
	
	design1 <- model.matrix(as.formula(opt$d1), data=mdata)
	
	if (opt$d1 != "~1") {
		design1 <- design1[match(colnames(m), rownames(design1)),]
	}
	
	v <- voom(M, design1, plot=FALSE)

	inputToBatchRm <- v
}



# **********************
#   removeBatchEffect
# **********************

design2 <- model.matrix(as.formula(opt$d2), data=mdata)

if (opt$d2 != "~1") {
	design2 <- design2[match(colnames(m), rownames(design2)),]
}

batch = mdata[match(colnames(m), mdata[, opt$merge_mdata_on]), opt$batch]

if (opt$method == "voom" | opt$method == "limma") {
	out = removeBatchEffect(inputToBatchRm, batch=batch, design=design2)
}

# Convert back to read counts
if (opt$read_counts) {
	out = pmax(sweep(2**out, 2, M$samples$norm.factors * M$samples$lib.size, FUN="*")/1e+06 - 0.5, 0)
}
	
	
# ****************
#    ComBat
# ****************     


if (opt$method == "combat") {

	suppressPackageStartupMessages(library(sva))
	out = ComBat(dat=m, batch=batch, mod=design2)

}


# =================== OUTPUT ======================

out = round(out, digits=5)
outF = ifelse(opt$output=="stdout", "", opt$output)
write.table(out, file=outF, quote=FALSE, sep="\t")

q(save='no')
