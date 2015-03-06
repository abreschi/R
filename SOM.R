#!/usr/bin/env Rscript 

options(stringsAsFactors=F)
set.seed(1)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin",
	help="Columns are samples and rows are dimensions. Can be stdin [default=%default]"),

make_option(c("-r", "--replace_NAs"), default=FALSE, action='store_true',
	help="Replace NAs by 0s. If FALSE, rows containing NAs are omitted. [default=%default]"),

make_option(c("-G", "--grid"), default="5x4",
	help="Grid structure. Format: gridRowsxgridCols. [default=%default]"),

make_option(c("-t", "--topology"), default="hexagonal",
	help="Grid topology. <hexagonal> | <rectangular> [default=%default]"),

make_option(c("-T", "--toroidal"), default=FALSE, action="store_true",
	help="Toroidal strucure [default=%default]"),

make_option(c("-n", "--iterations"), default=100, 
	help="Number of iteration [default=%default]"),

make_option(c("-m", "--metadata"), 
	help="tsv file with metadata on matrix experiment"),

make_option(c("-f", "--fields"), 
	help="choose the fields you want to use for super-organized maps, comma-separated. Needs a metadata file"),

make_option(c("-o", "--output"), default="SOM.out.tsv",
	help="Output file name, with extension .tsv [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="verbose output")
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	description = "
	Apply SOM to a matrix of values.
	The resulting cluster unit is appended as last column to the original matrix.",
	option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


##------------
## LIBRARIES
##------------ 

if(opt$verbose) {cat('Libraries loading... ')}

suppressPackageStartupMessages(library('reshape2'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library("kohonen"))

if (opt$verbose) {cat('DONE\n\n')}

##-------------##
## BEGIN       ##
##-------------##

# read the matrix from the command line
if (opt$input_matrix == "stdin") {
	input = "stdin"
} else {
	input = opt$input_matrix
}

#input = ifelse(opt$input_matrix == "stdin", file("stdin"), opt$input_matrix)

data = read.table(input, h=T)
init_nrow = nrow(data)
if (opt$verbose) {print(head(data))}

# Deal with NAs
if (opt$replace_NAs) {
	m = replace(data, is.na(data), 0)
} else {
	keep_index = which(rowSums(is.na(data)) == 0)
	m = na.omit(data)
	if(opt$verbose) {cat("Dimensions retained after omittig NAs:", nrow(m), "of", init_nrow, "\n")}
}

# Initialize the column with the clusters
data$K <- NA

# read the SOM grid 
gridStruct = as.numeric(strsplit(opt$grid, "x")[[1]])
gridRows = gridStruct[1]
gridCols = gridStruct[2]
if(opt$verbose) {cat("SOM grid:", gridStruct, "\n")}


# read the metadata from the metadata file if provided
if (!is.null(opt$metadata)) {
	merge_mdata_on = 'labExpId'
	toSuperSOM = list()
	mdata = read.table(opt$metadata, h=T, sep='\t')
	mdata[merge_mdata_on] <- gsub(",", ".", mdata[,merge_mdata_on])
	# Take only mdata rows which are in the columns of the input matrix
	mdata = mdata[which(mdata[,merge_mdata_on] %in% colnames(data)),]
	fields = strsplit(opt$fields, ",")[[1]]
	mdata = unique(mdata[,c(merge_mdata_on, fields)])
	levs = levels(as.factor(mdata[,fields]))
	for (lev in levs) {
		ids = mdata[which(mdata[, fields] == lev), merge_mdata_on]
		toSuperSOM[[lev]] <- as.matrix(m[ids])
	}
	SOM = supersom(
		toSuperSOM,
		grid = somgrid(gridCols, gridRows, opt$topology), 
		toroidal=opt$toroidal,
		rlen = opt$iterations
	)
	df_changes = melt(as.matrix(SOM$changes), varnames=c("x", "f"), value.name="y")
	data[keep_index, "K"] <- SOM$unit.classif
} else {
	SOM = som(
		as.matrix(m), 
		grid = somgrid(gridCols, gridRows, opt$topology), 
		toroidal = opt$toroidal,
		rlen = opt$iterations
	)
	df_changes = data.frame(x = seq_along(SOM$changes), y = SOM$changes)
	data[keep_index, "K"] <- SOM$unit.classif
}



# ~~~~~~~~~~~~ #
# OUTPUT       #
# ~~~~~~~~~~~~ #

write.table(data, file=opt$output, quote=FALSE, sep="\t")

# ~~~~~~~~~~ #
# PLOT       #
# ~~~~~~~~~~ #

out.prefix = gsub(".tsv", "", opt$output)

theme_set(theme_bw(base_size=16))

# convergence
gp = ggplot(df_changes, aes(x, y)) 
if (!is.null(opt$metadata)) {
	gp = gp + geom_line(aes(color=f, group=f))
} else {
	gp = gp + geom_line()
}
gp = gp + labs(x='Iteration', y='Mean distance to unit')
ggsave(sprintf('%s.convergence.pdf', out.prefix), h=5, w=5)

# Number of elements per unit
df = setNames(as.data.frame(table(SOM$unit.classif)), c('x', 'y'))
gp = ggplot(df, aes(x, y))
gp = gp + stat_identity(geom='bar', color='blue', fill='white')
gp = gp + labs(x='Unit', y='Number of elements')
gp = gp + theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(sprintf('%s.units.pdf', out.prefix), h=5, w=8)


q(save='no')
