#!/usr/bin/env Rscript 


# -- Variables --

options(stringsAsFactors=F)


##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-l", "--lcol"), help="Comma-separeted colors for the lines of the sets. Only names accepted for the moment! [default: black]"),
make_option(c("-f", "--fcol"), help="Comma-separated colors for the surfaces of the sets. Only names accepted for the moment! [default:palette]"),
make_option(c("-L", "--Lcol"), help="Comma-separated colors for the labels of the sets. Only names accepted for the moment! [default:black]"),

make_option(c("-W", "--width"), default=3, type='integer',
	help="width of the plot in inches [default=%default]"),

make_option(c("-H", "--height"), default=3, type='integer',
	help="height of the plot in inches [default=%default]"),

make_option(c("-o", "--output"), help="output file name WITHOUT extension [default=Venn.out]", default="venn.out")
)

cat("\nNOTE: NAs are treated as strings and not as missing values\n\n", file=stderr())

parser <- OptionParser(usage = "%prog [options] file(s)", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
arg <- arguments$args 
#print(arguments)



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
	venn_list = modifyList(venn_list, l)
}

for (i in seq_along(arg)) {
	l = as.list(read.table(arg[i], h=T, na.strings=NULL))
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


# change the label colors
if (is.null(opt$Lcol)) {
	label_col = rep('black', length(venn_list))
}else{
	label_col <- strsplit(opt$Lcol, ',')[[1]]  
}


# ===============
# ** plotting...
# ===============

image_type = 'png'
w = opt$width
h = opt$height


# Set the label distances
cat.dist = (w*c(0.01,0.01,0,0,0))[1:length(venn_list)]

if (length(venn_list) == 4) {
	cat.dist = (w*c(0.08, 0.08, 0.02, 0.02))
}

# Set the label position
if (length(venn_list) == 5) {
	cat.pos = c(20, 40, 40, 30, 10)
}

if (length(venn_list) == 4) {
	cat.pos = c(180, 180, 0, 0)
}

if (length(venn_list) == 3) {
	cat.pos = c(-20, 20, 20)
}
	
if (length(venn_list) == 2) {
	cat.pos = c(20, 30)
}



venn.diagram(venn_list, 
	filename=sprintf("%s.%s", opt$output, image_type), 
	imagetype=image_type,
	units='in',
	width=w,
	height=h,
	col=col,
	cat.col = label_col,
	cat.pos = cat.pos,
	cat.dist = cat.dist,
	cat.cex=c(0.8,0.8,0.8,0.8,0.8)[1:length(venn_list)],
	fill=face_col
)

# writing the intersection
write.table(data.frame(merged), file=sprintf("%s.tsv", opt$output), quote=F, row.names=F)

q(save='no')


