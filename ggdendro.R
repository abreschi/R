#!/usr/bin/env Rscript

# DEBUG OPTIONS
opt = list()
opt$input_matrix = "~/Documents/blueprint/pilot/Flux/Long/bp.human.long.gene.RPKM.idr_01.thr_0.names_False.tsv"
opt$col_metadata = "~/Documents/blueprint/pilot/bp_rna_dashboard_mmaps.crg.tsv"
opt$colSide_by = "cell"
opt$row_metadata = "/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/Long/gen15.gene.super.gene_type.with_header.tsv"
opt$merge_row_mdata_on = "gene"
opt$col_dendro = TRUE
opt$row_dendro = TRUE

# DEFAULT OPTIONS
opt$log = FALSE
opt$colSide_by = NULL
opt$col_labels = NULL
opt$row_labels = NULL
opt$merge_col_mdata_on = "labExpId"
opt$dist = "euclidean"
opt$hclust = "complete"
opt$base_size = 16


##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input_matrix"), 
	help="the matrix you want to analyze. \"stdin\" to read from standard input"),

make_option(c("-l", "--log"), action="store_true", default=FALSE, 
	help="apply the log10. NAs are treated as 0s and a pseudocount is added if specified [default=%default]"),

make_option(c("-p", "--pseudocount"), type="double", default=1e-04,
	help="specify a pseudocount for the log [default=%default]"),

make_option(c("--col_metadata"), 
	help="one tsv file with metadata on matrix columns. Can be left empty."),

make_option(c("--merge_col_mdata_on"), default="labExpId",
	help="which field of the metadata corresponds to the column headers? [default=%default]"), 

#make_option(c("--row_metadata"), 
#	help="one tsv file with metadata on matrix rows. Can be left empty."),
#
#make_option(c("--merge_row_mdata_on"), default="labExpId",
#	help="which field of the metadata corresponds to the row names? [default=%default]"), 
#
make_option(c("--col_labels"), 
	help="Specify the field for the col labels. \"none\" for no col labels. If empty the column headers are used."),

#make_option(c("--row_labels"), 
#	help="Specify the field for the col labels. \"none\" for no row labels. If empty the row names are used."),
#
make_option(c("--colSide_by"), 
	help="Specify the field(s), you want the column sides coloured by. If empty no color side is added."),

#make_option(c("--rowSide_by"), 
#	help="Specify the field(s), you want the row sides coloured by. If empty no color side is added."),
#
#make_option(c("--col_dendro"), action="store_true", default=FALSE, 
#	help="Print the column dendrogram [default=%default]"),
#
#make_option(c("--row_dendro"), action="store_true", default=FALSE, 
#	help="Print the row dendrogram [default=%default]"),
#
make_option(c("-V", "--vertical"), action="store_true", default=FALSE,
	help="Draw the cluster vartically [default=%default]"),

make_option(c("-d", "--dist"), default="euclidean",
	help="distance measure between columns. Choose among <p> (pearson), <s> (spearman),
		all methods supported by the function dist(). [default=%default]"),

make_option(c("-c", "--hclust"), default="complete",
	help="hierarchical clustering method. Choose among the method of the function hclust(). [default=%default]"),

make_option(c("-B", "--base_size"), default=16,
	help="The font base size as defined in ggplot2. [default=%default]"),

make_option(c("-W", "--width"), type="integer",
	help="Choose the heatmap width in inches. Default is proportional to the number of columns"),

make_option(c("-H", "--height"), type="integer",
	help="Choose the heatmap height in inches. Default is proportional to the number of rows"),

make_option(c("-o", "--output"), 
	help="Output file name, with the extension. [default=%default]", default="ggdendro.out.pdf")

)


parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
#print(opt)


#------------#
# LIBRARIES  #
#------------#

cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(grid))
cat("DONE\n\n")


# ==========================================
# Function for extracting legend from ggplot
# ==========================================

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}


# ======================
# Plotting variables
# ======================

base_size = opt$base_size
theme_set(theme_grey(base_size))
theme_update(axis.ticks=element_blank())
theme_update(axis.ticks.margin = unit(0.01, "inch"))
theme_update(axis.ticks.length = unit(0.01, "inch"))



# ===== #
# BEGIN #
# ===== #


# read table
if (opt$input_matrix == "stdin") {
	m = read.table(file("stdin"), h=T)} else {
	m = read.table(opt$input_matrix, h=T)
}


m = m[1:1000,]

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# apply the log10 if needed
if (opt$log) {m <- log10(replace(m, is.na(m), 0) + opt$pseudocount)}

# melt the data frame
df = melt(as.matrix(m))


# --------------- Metadata processing -------------

# read metadata
if (!is.null(opt$col_metadata)) {col_mdata = read.table(opt$col_metadata, h=T, sep="\t")}
#if (!is.null(opt$row_metadata)) {row_mdata = read.table(opt$row_metadata, h=T, sep="\t")}
# read which fields are needed from the metadata
if (!is.null(opt$colSide_by)) {colSide_by = strsplit(opt$colSide_by, ",")[[1]]} else {colSide_by = NULL}
#if (!is.null(opt$rowSide_by)) {rowSide_by = strsplit(opt$rowSide_by, ",")[[1]]} else {rowSide_by = NULL}
if (!is.null(opt$col_labels) && opt$col_labels != "none") {col_label_fields = strsplit(opt$col_labels,",")[[1]]} else {col_label_fields=NULL}
#if (!is.null(opt$row_labels) && opt$row_labels != "none") {row_label_fields = strsplit(opt$row_labels,",")[[1]]} else {row_label_fields=NULL}

col_mdata_header = unique(c(opt$merge_col_mdata_on, colSide_by, col_label_fields))
#row_mdata_header = unique(c(opt$merge_row_mdata_on, rowSide_by, row_label_fields))

# merge metadata and data (NB: The column Var2 stays)
if (!is.null(opt$col_metadata)) {
	col_mdata[opt$merge_col_mdata_on] <- gsub(",", ".", col_mdata[,opt$merge_col_mdata_on])
	df = merge(df, col_mdata[col_mdata_header], by.x="Var2", by.y=opt$merge_col_mdata_on)
}
#if (!is.null(opt$row_metadata)) {
#	row_mdata[opt$merge_row_mdata_on] <- gsub(",", ".", row_mdata[,opt$merge_row_mdata_on])
#	df = merge(df, row_mdata[row_mdata_header], by.x="Var1", by.y=opt$merge_row_mdata_on)
#}


print(head(df))

# ---------------- Dendrogram ----------------------

# COLUMNS

if (opt$vertical) {
	dendro_hjust=0
} else {
	dendro_hjust=1
	dendro_angle=90
}

if (opt$dist == "p" || opt$dist =="s") {
	colDist = as.dist(1-cor(m, method=opt$dist, use="p"))
} else {
	colDist = dist(t(m), method=opt$dist)
}
colHC = hclust(colDist, method=opt$hclust)
colHC_data = dendro_data(as.dendrogram(colHC))
col_ggdendro = ggplot(segment(colHC_data))
col_ggdendro = col_ggdendro + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
#col_ggdendro = col_ggdendro + geom_text(data=label(colHC_data), 
#	aes(x=x, y=y, label=label), angle=dendro_angle, hjust=dendro_hjust) 
if (opt$vertical) {
	col_ggdendro = col_ggdendro + coord_flip()
} 
col_ggdendro = col_ggdendro + scale_x_continuous(expand=c(0.0, 0.5), labels=NULL) 
col_ggdendro = col_ggdendro + scale_y_continuous(expand=c(0.0, 0.0), labels=NULL)
col_ggdendro = col_ggdendro + theme(plot.margin=unit(c(0.10, 0.00, 0.00, 0.01), "inch")) # top, right, bottom, left
col_ggdendro = col_ggdendro + theme_dendro()
col_ggdendro = col_ggdendro + labs(x=NULL, y=NULL)
#col_ggdendro


# -------------------- Column Side Colors ------------

ColSides = list(); ColSide_legends = list()
col_limits = label(colHC_data)[,3]

print(col_limits)
if (!is.null(opt$colSide_by)) {
	i=1;
	for (colSide in colSide_by) {
		colSide_data = unique(df[c("Var2", colSide)])
		ColSide = ggplot(colSide_data, aes(x=Var2, y="a"))
		ColSide = ColSide + geom_tile(aes_string(fill=colSide), color="black")
		ColSide = ColSide + scale_x_discrete(limits = col_limits, labels=NULL, expand=c(0,0))
		ColSide = ColSide + scale_y_discrete(labels=NULL, expand=c(0,0))
		ColSide = ColSide + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.01),"inch"))
		ColSide = ColSide + labs(x=NULL, y=NULL)
		ColSide_legends[[i]] = g_legend(ColSide)
		ColSide = ColSide + theme(legend.position="none")
		ColSides[[i]] = ColSide; i=i+1;
	}
}


#ggsave('test.pdf', plot=ColSides[[1]], h=10, w=7)
#q(save='no')


# This works in the X11 device
#row_labels_inches = max(strwidth(row_labels, units="in"))
#col_labels_inches = max(strwidth(col_labels, units="in"))

# This works in the pdf device
#row_labels_inches = max(strwidth(row_labels, units="in", cex=base_size*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))
#col_labels_inches = max(strwidth(col_labels, units="in", cex=base_size*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))





# ============================
# Compose with viewports
# ============================


# >>> Column dendrogram viewport <<<<<<<<<

dendro_vp_y = 0.25*length(ColSides) # + col_labels_inches
dendro_vp_x = 0.1
dendro_h = 3
dendro_w = 3.5 

colDendro_vp = viewport(
y = dendro_vp_y, 
x = dendro_vp_x,
h = dendro_h,
w = dendro_w,
default.units = "inch",
just = c("left", "bottom")
)


#total_h = matrix_vp_y + matrix_vp_h + 0.25*length(colSide_by) + max(dendro_h, matrix_scale_h)
#total_w = matrix_vp_x + matrix_vp_w + 0.25*length(rowSide_by) + max(row_dendro_w, matrix_scale_w)

total_h = dendro_h + dendro_vp_y
total_w = dendro_w



# >>>>>> Column side viewport <<<<<<<<<<<

if (!is.null(opt$colSide_by)){
	ColSide_vps = list(); ColSide_label_vps = list()
	for (i in 1:length(ColSides)) {
		ColSide_vps[[i]] = viewport(
			y = 0.01 + 0.25*(i-1), 
			x = dendro_vp_x,
			h = 0.25,
			w = dendro_w,
			default.units = "inch",
			just = c("left", "bottom") 
		 )
		ColSide_label_vps[[i]] = viewport(
			y = 0.01 + 0.25*(i-1),
			x = dendro_vp_x + dendro_w,
			h = 0.25,
			w = as.numeric(strwidth(colSide_by, "inch")),
			default.units = "inch",
			just = c("left", "bottom")
		)
	}
	total_w = total_w + 1.5 
}





# >>>>> Column and row side legends viewport <<<<<<<<<<<<<<<<<


if (!is.null(opt$colSide_by)) {
	legend_width_inch = max(sapply(ColSide_legends, function(x) sum(sapply(x$widths, convertUnit, "in"))))
	legend_height_inch = sum(sapply(ColSide_legends, function(x) sum(sapply(x$heights, convertUnit, "in")))) + 0.1*(length(ColSide_legends)-1)
#	legend_height_inch = max(sapply(ColSide_legends, function(x) sum(sapply(x$heights, convertUnit, "in"))))
#	legend_width_inch = max(strwidth(unlist(unique(df[c(colSide_by, rowSide_by)])), unit="inch")) + 0.4
#	legend_height_inch = total_h/length(c(colSide_by, rowSide_by))
	side_legend_vps = list()
	for (i in 1:length(colSide_by)) {
		side_legend_vp = viewport(
			y = dendro_vp_y + dendro_h + 0.25*length(colSide_by) - legend_height_inch*(i-1),
			x = dendro_vp_x + dendro_w,
			h = legend_height_inch,
			w = legend_width_inch,
			default.units = "inch",
			just = c("left", "top")
		)
		side_legend_vps[[i]] = side_legend_vp
	}
} else { legend_width_inch = 0}

#total_w = total_w + legend_width_inch
#if (!is.null(opt$colSide_by) || !is.null(opt$rowSide_by)) {
#	legend_width_inch = max(strwidth(unlist(unique(df[c(colSide_by, rowSide_by)])), unit="inch")) + 0.4
#	legend_height_inch = total_h/length(c(colSide_by, rowSide_by))
#	side_legend_vps = list()
#	for (i in 1:length(c(colSide_by, rowSide_by))) {
#		side_legend_vp = viewport(
#			y = matrix_vp_y + matrix_vp_h + 0.25*length(colSide_by) + col_dendro_h - legend_height_inch*(i-1),
#			x = matrix_vp_x + matrix_vp_w + 0.25*length(rowSide_by) + row_dendro_w,
#			h = legend_height_inch,
#			w = legend_width_inch,
#			default.units = "inch",
#			just = c("left", "top")
#		)
#		side_legend_vps[[i]] = side_legend_vp
#	}
#} else { legend_width_inch = 0}
#
#total_w = total_w + legend_width_inch


# =======================================
# PRINT PLOT
# =======================================

pdf(opt$output, h = total_h, w=total_w)

#X11(h=total_h, w=total_w)

## Print matrix
#print(p1, matrix_vp, newpage=FALSE)
#
## Print matrix legend
#pushViewport(matrix_scale_vp); grid.draw(p1_legend); upViewport()

# Print column side colors
if (!is.null(opt$colSide_by)) {
	for (i in 1:length(ColSide_vps)) {
		print(ColSides[[i]], vp=ColSide_vps[[i]], newpage=FALSE)
	#	grid.text(colSide_by[i], x=unit(1,"npc"), just="left", vp=ColSide_label_vps[[i]], gp=gpar(face="bold"))
	}
}

## Print row side colors
#if (!is.null(opt$rowSide_by)) {
#	for (i in 1:length(RowSide_vps)) {
#		print(RowSides[[i]], vp=RowSide_vps[[i]], newpage=FALSE)
#		grid.text(rowSide_by[i], y=unit(1,"npc"), just="right", rot=90, vp=RowSide_label_vps[[i]], gp=gpar(face="bold"))
#	}
#}
#
## Print column dendrogram
#if (opt$col_dendro) {print(col_ggdendro, vp=colDendro_vp, newpage=FALSE)}
print(col_ggdendro, vp=colDendro_vp, newpage=FALSE)
#
## Print row dendrogram
#if (opt$row_dendro) {print(row_ggdendro, vp=rowDendro_vp, newpage=FALSE)}
#
# Print column side color scales
if (!is.null(opt$colSide_by)) {
	for (i in 1:length(colSide_by)) {
		all_side_legends = ColSide_legends
		pushViewport(side_legend_vps[[i]]); grid.draw(all_side_legends[[i]]); upViewport()
	}
}

## Print column and row side color scales
#if (!is.null(opt$colSide_by) || !is.null(opt$rowSide_by)) {
#	for (i in 1:length(c(colSide_by, rowSide_by))) {
#		all_side_legends = c(ColSide_legends, RowSide_legends)
#		pushViewport(side_legend_vps[[i]]); grid.draw(all_side_legends[[i]]); upViewport()
#	}
#}
#

dev.off()

file.remove("Rplots.pdf")
q(save="no")
