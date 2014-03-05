#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(igraph))

opt = list()
opt$input = "top500.sd.RPKM.glasso.tsv"
opt$node_size = 2

options(stringsAsFactors=FALSE)

set.seed(123)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin. Columns are node1, node2, weigth [default=%default]"),

make_option(c("-o", "--output"), default="network.pdf",
	help="Output file name [default=%default]"),

make_option(c("--nodes"), 
	help="File with node attributes"),

make_option(c("--node_color"), type="integer",
	help="Index of the node color"),

make_option(c("--label_color"), type="integer",
	help="Index of the label color"),

make_option(c("--node_frame_color"), type="integer",
	help="Index of the node frame color"),

make_option(c("--node_size"), default=2,
	help="Size of the nodes [default=%default]"),

make_option(c("--palette"), default="/users/rg/abreschi/R/palettes/rainbow.15.txt",
	help="File with colorname in RGB format [default=%default]"),

make_option(c("--label_cex"), default=1, type="double",
	help="Size of the labels in cex [default=%default]"),

make_option(c("--normalize"), action="store_true", default=FALSE,
	help="Normalize to have all 1s in the diagonal [default=%default]"),

make_option(c("--diag"), action="store_true", default=FALSE,
	help="Report also edges to and from the same node [default=%default]"),

make_option(c("--directed"), action="store_true", default=FALSE,
	help="Is the graph directed? [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


# BEGIN

# Read input
if (opt$input == "stdin") {
	m = read.table(file("stdin"), h=F)
} else {
	m = read.table(opt$input, h=F)
}

palette = read.table(opt$palette, h=F, comment.char="%")$V1


# Creat graph
g = graph.data.frame(m, directed=opt$directed)


# Normalize if needed
M = as.matrix(get.adjacency(g))
if (opt$normalize) {
	lambda = 1/sqrt(diag(M))
	M = sweep(sweep(M, MARGIN=2, lambda, `*`), MARGIN=1, lambda, `*`)
	g = graph.adjacency(as.matrix(M), mode=c("max"), weighted=TRUE, diag=opt$diag)
}


# Get node attributes

V(g)$color = rep("white", vcount(g))
V(g)$frame.color = rep("white", vcount(g))
V(g)$label.color = rep("black", vcount(g))
V(g)$size = rep(opt$node_size, vcount(g))
V(g)$label.cex = rep(opt$label_cex, vcount(g))


if (!is.null(opt$nodes)) {
	node_attr = read.table(opt$nodes, h=F)
	match_node = match(node_attr[,1], V(g)$name)

	# Get node color
	if (!is.null(opt$node_color)) {
		node_colors = palette[node_attr[,opt$node_color]]
		V(g)[match_node]$color = node_colors
	}
	
	# Get label color
	if (!is.null(opt$label_color)) {
		label_colors = palette[node_attr[,opt$label_color]]
		V(g)[match_node]$label.color = label_colors
	}

	# Get node frame color
	if (!is.null(opt$node_frame_color)) {
		node_frame_colors = palette[node_attr[,opt$node_frame_color]]
		V(g)[match_node]$frame.color = node_frame_colors
	}
}



# PLOT

pdf(opt$output, h=7, w=7)

plot(
	g, 
	layout=layout.fruchterman.reingold, 
#	vertex.label=NA, 
	main=""
)

if (!is.null(opt$nodes) & !is.null(opt$node_color)) {
	labels = unique(node_attr[,opt$node_color])
	legend(
		"topright", 
		legend=labels, 
		fill=node_colors[match(labels, node_attr[,opt$node_color])]
	)
}

dev.off()


q(save='no')

