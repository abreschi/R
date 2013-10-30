##------------
## LIBRARIES
##------------

cat("Loading libraries... ") 
#suppressPackageStartupMessages(library(reshape2))
#suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
#suppressPackageStartupMessages(library(plyr))
cat("DONE\n\n")

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-A", "--matrix_A"), help="the matrix you want to subtract from, WITH header (A-B)"),
make_option(c("-B", "--matrix_B"), help="the matrix you want to subtract, WITH header (A-B)"),
make_option(c("-r", "--replace_NA_with"), help="value you want to replace NA with, if null NAs are not replaced and difference will be NA", type="numeric"),
make_option(c("-o", "--output_prefix"), help="additional prefix for otuput [default=%default]", default='out'),
make_option(c("-e", "--expression"), help="expression you want to corresponding cells of the matrices, e.g. A-B")

#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
#make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
#make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
#make_option(c("-f", "--fill_by"), help="choose the color you want to fill by [default=NA]", type='character', default=NA),
#make_option(c("-t", "--tags"), help="comma-separated field names you want to display in the labels", default="cell,sex,age")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)

###################
# BEGIN           #
###################

A = read.table(opt$matrix_A, h=T)
B = read.table(opt$matrix_B, h=T)

if (!is.null(opt$replace_NA_with)) {
A <- replace(A, is.na(A), opt$replace_NA_with)
B <- replace(B, is.na(B), opt$replace_NA_with)
}

M = eval(parse(text=opt$expression))

output = sprintf("%s.tsv", opt$output_prefix)
write.table(M, output, quote=F, row.names=T, sep="\t")

q(save='no')
