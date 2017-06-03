#!/usr/bin/env Rscript

set.seed(1555)

# DEFAULT OPTIONS

opt = list()
opt$log10 = FALSE
opt$pseudocount = 1e-04
opt$row_as_variables = FALSE

##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin", help="Matrix of observations (rows are the obs) or distances [default=%default]"),
make_option(c("-l", "--log10"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-p", "--pseudocount"), type="double", help="specify a pseudocount for the log [default=%default]", default=0),

make_option(c("-T", "--input_type"), default="dist", help="Choose between dist and obs. [default=%default]"),

make_option(c("-n", "--iter"), default=1000, help="Number of iterations [default=%default]"),

make_option(c("-m", "--metadata"), help="Tab-separated file with the metadata. Can be omitted"),

make_option(c("-M", "--method"), help="Choose one of: <MDS> <tSNE>"),

make_option(c("--merge_mdata_on"), default="labExpId",
	help="[default=%default]"),

make_option(c("--perplexity"), default=30, help="[default=%default]"),

#make_option(c("-o", "--output"), help="additional info you want to put in the output file name", default="out"),
make_option(c("-c", "--color_by"), help="choose the fields in the metadata you want to color by", type='character'),

make_option(c("-s", "--shape_by"), help="choose the fields in the metadata you want to shape by", type='character'),

make_option(c("-P", "--palette"), help="palette file [default=%default]"),

make_option(c("--shapes"), 
	help="File with the shapes [default=%default]"),

make_option(c("-b", "--base-size"), dest="base_size", default=18,
	help="Text size [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help="Width of the plot in inches [default=%default]"),

make_option(c("-H", "--height"), default=7,
	help="Height of the plot in inches [default=%default]"),

make_option(c("-o", "--output"), default="Rtsne",
	help="Prefix for the file names. [default=%default]"),

make_option(c("-v", "--verbose"), action='store_true', default=FALSE,
	help="verbose output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list,
	description="\nPerform Barnes Hut t-SNE. Return a plot with the points in the new 2D and a plot of the cost")
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if (opt$verbose) {print(opt)}

suppressPackageStartupMessages(library("reshape"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("Rtsne"))
suppressPackageStartupMessages(library("plyr"))

###############
# BEGIN
##############

# read input matrix
if (opt$input_matrix == "stdin") {input=file("stdin")} else {input=opt$input_matrix}
m = read.table(input, h=T, sep="\t")
id = colnames(m)

# --- dist ---
if (opt$input_type == "dist") {
	m = as.dist(m)
}

# --- obs ---
if (opt$input_type == "obs") {
	# apply the log if required
	if (opt$log10) {m = log10(replace(m, is.na(m), 0) + opt$pseudocount)}
	m = t(m)
}


## Read the color palette
if (!is.null(opt$palette)) {my_palette = read.table(opt$palette, h=F, comment.char="%", sep="\t")$V1}


# Read the shapes
if (!is.null(opt$shapes)) {
	my_shapes = read.table(opt$shapes, h=F, comment.char="%")$V1
}

if (opt$method == "tSNE") {
	tsne = Rtsne(m, max_iter=opt$iter, perplexity=opt$perplexity, pca=FALSE, is.distance=(opt$input_type=="dist"))
	x = tsne$Y[,1]
	y = tsne$Y[,2]

	# Plot the cost at each iteration
	cost_df = data.frame(iter=seq_along(tsne$itercosts), cost=tsne$itercosts)
	gp = ggplot(cost_df, aes(x=iter, y=cost)) + geom_line()
	ggsave(sprintf("%s.cost.pdf", opt$output), h=5, w=7)

}

if (opt$method == "MDS") {
	fit <- cmdscale(m, eig=TRUE, k=2) 
	x <- fit$points[,1]
	y <- fit$points[,2]
}


# ***************
# Plot
# ***************

map_df = data.frame(x=x, y=y, id=id)


# HANDLING METADATA

if (!is.null(opt$metadata)){
	mdata = read.table(opt$metadata, h=T, sep="\t", row.names=NULL, quote="", check.names=F);
	mdata[,opt$merge_mdata_on] <- gsub("[,-]", ".", mdata[,opt$merge_mdata_on])
	mdata <- unique(mdata[,c(opt$merge_mdata_on, opt$color_by, opt$shape_by)])
	map_df <- merge(map_df, mdata, by.x="id", by.y=opt$merge_mdata_on)
}


# PLOTTING

theme_set(theme_bw(base_size=opt$base_size))
theme_update(
	legend.title = element_text(size=opt$base_size),
	legend.key = element_blank(),
	panel.grid = element_blank()
)


geom_params = list()
mapping = list()

geom_params$size = 4

if (is.null(opt$color_by)) {
	geom_params$color = "black"
} else {
	mapping = modifyList(mapping, aes_string(color=paste("`", opt$color_by, "`", sep="")))
}

if (is.null(opt$shape_by)) {
	geom_params$shape = 19
} else {
	mapping = modifyList(mapping, aes_string(shape=paste("`", opt$shape_by, "`", sep="")))
}

class(mapping) <- "uneval"

pointLayer <- layer(
	geom = "point",
#	geom_params = geom_params,
	params = geom_params,
	stat = "identity",
	position = "identity",
	mapping = mapping	
)


gp = ggplot(map_df, aes(x=x, y=y)) + pointLayer
gp = gp + labs(x=paste(opt$method, "1", sep=" "), y=paste(opt$method, "2", sep=" "), title="")
if (!is.null(opt$palette)) {
	if (is.discrete(map_df[,opt$color_by])) {
		gp = gp + scale_color_manual(values=my_palette)
	} else {
		gp = gp + scale_color_gradientn(colours=my_palette)
	}
}
if (!is.null(opt$shapes)) {
	gp = gp + scale_shape_manual(name=opt$shape_by, values=my_shapes);
}


ggsave(sprintf("%s.map.pdf", opt$output), w=opt$width, h=opt$height)

q(save='no')
