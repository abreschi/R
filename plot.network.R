#!/usr/bin/env Rscript


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
	help="File with node attributes, no header"),

make_option(c("--node_color"), type="integer",
	help="Index of the node color"),

make_option(c("--label"),
	help="Leave empty for using the node names. \"none\" for no label."),

make_option(c("--label_color"), type="integer",
	help="Index of the label color"),

make_option(c("--node_frame_color"), type="integer",
	help="Index of the node frame color"),

make_option(c("--node_size"), default=2,
	help="Size of the nodes [default=%default]"),

make_option(c("--node_shape"), type='integer',
	help="Index of the node shape"),

make_option(c("--node_palette"), default="/users/rg/abreschi/R/palettes/rainbow.15.txt",
	help="File with colorname in RGB format [default=%default]"),

make_option(c("--node_frame_palette"), default="/users/rg/abreschi/R/palettes/rainbow.4.txt",
	help="File with colorname in RGB format [default=%default]"),

make_option(c("--label_cex"), default=1, type="double",
	help="Size of the labels in cex [default=%default]"),

make_option(c("--normalize"), action="store_true", default=FALSE,
	help="Normalize to have all 1s in the diagonal [default=%default]"),

make_option(c("--diag"), action="store_true", default=FALSE,
	help="Report also edges to and from the same node [default=%default]"),

make_option(c("--directed"), action="store_true", default=FALSE,
	help="Is the graph directed? [default=%default]"),

make_option(c("-H", "--height"), default=9,
	help="Height of the plot in inches [default=%default]"),

make_option(c("-W", "--width"), default=9,
	help="Width of the plot in inches [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(igraph))


#################################################################
# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
 
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add.vertex.shape("triangle", clip=vertex.shapes("circle")$clip,
                 plot=mytriangle)



# BEGIN

# Read input
if (opt$input == "stdin") {
	m = read.table(file("stdin"), h=F)
} else {
	m = read.table(opt$input, h=F)
}

node_palette = read.table(opt$node_palette, h=F, comment.char="%")$V1
node_frame_palette = read.table(opt$node_frame_palette, h=F, comment.char="%")$V1
node_shape_palette = c("circle", "square", "csquare", "rectangle", "crectangle", "vrectangle", "none")

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
V(g)$shape = rep("circle", vcount(g))
V(g)$size = rep(opt$node_size, vcount(g))
V(g)$label.cex = rep(opt$label_cex, vcount(g))


if (!is.null(opt$nodes)) {
	node_attr = read.table(opt$nodes, h=F, quote="\"")
	# Select only metadata rows for nodes in the network
	node_attr = node_attr[node_attr[,1] %in% V(g)$name,]
	match_node = match(V(g)$name, node_attr[,1])

	# Get node color
	if (!is.null(opt$node_color)) {
		node_colors = node_palette[as.factor(node_attr[,opt$node_color])]
		V(g)$color = node_colors[match_node]
	}

	# Get label color
	if (!is.null(opt$label_color)) {
		label_colors = palette[node_attr[,opt$label_color]]
		V(g)$label.color = label_colors[match_node]
	}

	# Get node frame color
#	if (!is.null(opt$node_frame_color)) {
#		node_frame_colors = node_frame_palette[as.factor(node_attr[,opt$node_frame_color])]
#		V(g)$frame.color = node_frame_colors[match_node]
#	}

	# Get label
	if (!is.null(opt$label) && opt$label != "none") {
		opt$label = as.integer(opt$label)
		node_labels = node_attr[,opt$label]
		V(g)$label = node_labels[match_node]
	}
	
	# Get node shape
	if (!is.null(opt$node_shape)) {
		node_shapes = node_shape_palette[as.factor(node_attr[,opt$node_shape])]
		V(g)$shape = node_shapes[match_node]
	}

}


vertex.label=NULL
if (opt$label == "none") {V(g)$label = NA}

#E(g)$label = m[,4]

# PLOT

pdf(opt$output, h=opt$height, w=opt$width)

cat("PLOTTING...")

plot(
	g, 
	layout=layout.kamada.kawai, 
#	layout=layout.reingold.tilford(g, root="1")
#	vertex.label=vertex.label, 
)

# Node color legend
if (!is.null(opt$nodes) & !is.null(opt$node_color)) {
	labels = unique(node_attr[,opt$node_color])
	legend(
		"topright", 
		legend=labels, 
		fill=node_colors[match(labels, node_attr[,opt$node_color])]
	)
}


# Node frame color legend
if (!is.null(opt$nodes) & !is.null(opt$node_frame_color)) {
	labels = unique(node_attr[,opt$node_frame_color])
	legend(
		"topleft", 
		legend=labels, 
		fill="white",
		border=node_frame_colors[match(labels, node_attr[,opt$node_frame_color])]
	)
}

cat(" DONE\n")

dev.off()


q(save='no')

