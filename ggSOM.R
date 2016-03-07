#!/usr/bin/env Rscript 

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input"), default="stdin",
	help="Tab-delimited file with unit id and attribute, no header. Can be stdin [default=%default]"),

#make_option(c("-r", "--replace_NAs"), default=FALSE, action='store_true',
#        help="Replace NAs by 0s. If FALSE, rows containing NAs are omitted. [default=%default]"),
#
make_option(c("-G", "--grid"), default="5x4",
	help="Grid structure. Format: gridRowsxgridCols. [default=%default]"),

make_option(c("-k", "--cluster"), type='integer',
	help="Index of the column with the cluster id"),

make_option(c("-f", "--factor"), type='integer',
	help="Index of the column with the factor"),

make_option(c("--func"), default="median",
	help="Function to aggregate in each cluster [default=%default]"),

make_option(c("-d", "--is_discrete"), default=FALSE, action="store_true",
	help="The factor is a discrete class. [default=%default]"),

make_option(c("-p", "--palette"), 
	help="File with custom palette. If no palette is given, the ggplot default is used."),

make_option(c("-T", "--title"),
	help="Main title"),

#make_option(c("-t", "--topology"), default="hexagonal",
#        help="Grid topology. <hexagonal> | <rectangular> [default=%default]"),
#
#make_option(c("-m", "--metadata"), 
#        help="tsv file with metadata on matrix experiment"),
#
#make_option(c("-f", "--fields"), 
#        help="choose the fields you want to use for super-organized maps, comma-separated. Needs a metadata file"),
#
make_option(c("-o", "--output"), default="ggSOM.pdf",
        help="Output file name [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE, 
	help="verbose output")
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	description = "
	Plot the output of SOM as hexagonal grid.
	",
	option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


##################
# LIBRARIES      #
##################

suppressPackageStartupMessages(library('ggplot2'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to compute the hexagon borders
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hex.x = function(lx) {
	c(
	lx + 0, 
	lx + 0.5, 
	lx + 0.5, 
	lx + 0, 
	lx - 0.5, 
	lx - 0.5 
	)
}

hex.y = function(ly) {
	c(
	ly - 1/(sqrt(3)), 
	ly - 1/(2*sqrt(3)), 
	ly + 1/(2*sqrt(3)), 
	ly + 1/(sqrt(3)), 
	ly + 1/(2*sqrt(3)), 
	ly - 1/(2*sqrt(3))
	)
}


##################
# BEGIN
##################


# read the SOM grid 
gridStruct = as.numeric(strsplit(opt$grid, "x")[[1]])
gridRows = gridStruct[1]
gridCols = gridStruct[2]
if(opt$verbose) {cat("SOM grid:", gridStruct, "\n")}

# Read the id attributes
if (opt$input == "stdin") {input = file('stdin')} else {input = opt$input}
m = read.table(input, h=F)

# Read the palette if provided
if (!is.null(opt$palette)) {
	palette = read.table(opt$palette, h=F, comment.char="%")[,1]
}

# Name the cluster column
colnames(m)[opt$cluster] <- 'id'

# Build the grid
hex = data.frame()
k=0

for (i in seq(gridRows)) {
	for (j in seq(gridCols)) {
		if (i%%2 != 0) {j = j+0.5}
		k = k+1
		hex = rbind(hex, data.frame(x = hex.x(j), y = hex.y(i), id=k))
	}
}

# Aggregate the attribute by cluster
factor_col = colnames(m)[opt$factor]


if (!opt$is_discrete) {
	formula_agg = as.formula(sprintf("%s~id", factor_col))
	m_agg = aggregate(formula_agg, m, eval(opt$func), na.rm=TRUE)
	# Merge the attribute with the hex grid
	df = merge(hex, m_agg)
} else {
	df = merge(hex, m)
}



########
# PLOT
########

#gp = ggplot(hex, aes(x, y)) + geom_polygon(aes(group=id, fill=id))
gp = ggplot(df, aes(x, y)) + geom_polygon(aes_string(group="id", fill=factor_col), color='black')
gp = gp + labs(title=opt$title, x=NULL, y=NULL)
if (!is.null(opt$palette)) {
	if (opt$is_discrete) {
		gp = gp + scale_fill_manual(values=palette)
	} else {
		gp = gp + scale_fill_gradientn(colours=palette)
	}
}

ggsave(opt$output, h=5, w=8)


