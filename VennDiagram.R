#!/usr/bin/env Rscript 


# -- Variables --

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-l", "--lcol"), help="Comma-separeted colors for the lines of the sets. Only names accepted for the moment! [default: black]"),
make_option(c("-f", "--fcol"), help="Comma-separated colors for the surfaces of the sets. Only names accepted for the moment! [default:palette]"),
make_option(c("-L", "--Lcol"), help="Comma-separated colors for the labels of the sets. Only names accepted for the moment! [default:black]"),
make_option(c("-o", "--output"), help="output file name WITHOUT extension [default=Venn.out]", default="venn.out")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
arg <- arguments$args 
print(arguments)



##------------
## LIBRARIES
##------------ 

cat("Loading libraries... ")
suppressPackageStartupMessages(library('VennDiagram'))
cat("DONE\n\n")


################
# BEGIN
################

# read the lists of elements from args
venn_list=list()
for (f in arg) {
l = as.list(read.table(f, h=T))
venn_list = c(venn_list, l)
}

for (i in seq_along(arg)) {
	l = as.list(read.table(arg[i], h=T))
	if (i==1) {
		venn_list = l
		merged = l[[1]]
		}
	if (i!=1) {
		venn_list = c(venn_list, l)
		merged = intersect(merged, l[[1]])
	}
}



# graphical parameters
#-----------------------

# change the line colors
if (is.null(opt$lcol)) {
	col <- rep('black', length(venn_list))
	
}else{
	col <- strsplit(opt$lcol, ',')[[1]]  
}


# change the surface colors
if (is.null(opt$fcol)) {
	face_col <- rainbow(length(venn_list))
}else{
	face_col = strsplit(opt$fcol, ',')[[1]]
}


#if (!is.null(opt$fcol)) {
#	face_col = strsplit(opt$fcol, ',')[[1]]
#}

# change the label colors
if (is.null(opt$Lcol)) {
	label_col = rep('black', length(venn_list))
}else{
	label_col <- strsplit(opt$Lcol, ',')[[1]]  
}



# plotting...
venn.diagram(venn_list, 
	filename=sprintf("%s.tiff", opt$output), 
	col=col,
	cat.col=label_col,
	fill=face_col
)

# writing the intersection
write.table(data.frame(merged), file=sprintf("%s.tsv", opt$output), quote=F, row.names=F)

q(save='no')


