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

make_option(c("-l", "--log"), 
	help="If you want to apply the log specify a base. Leave empty otherwise"),

make_option(c("-p", "--pseudocount"), type="double", default=1e-04,
	help="specify a pseudocount for the log [default=%default]"),

make_option(c("--col_metadata"), 
	help="one tsv file with metadata on matrix columns. Can be left empty."),

make_option(c("--merge_col_mdata_on"), default="labExpId",
	help="which field of the metadata corresponds to the column headers? [default=%default]"), 

make_option(c("--col_labels"), 
	help="Specify the field for the col labels. \"none\" for no col labels. If empty the column headers are used."),

make_option(c("--colSide_by"), 
	help="Specify the field(s), you want the column sides coloured by. If empty no color side is added."),

make_option(c("-P", "--colSide_palette"), 
	help="A palette file to color the side. Default is hue"),

make_option(c("-V", "--vertical"), action="store_true", default=FALSE,
	help="Draw the cluster vartically [default=%default]"),

make_option(c("-d", "--dist"), default="euclidean",
	help="distance measure between columns. Choose among <p> (pearson), <s> (spearman),
		all methods supported by the function dist(). [default=%default]"),

make_option(c("-c", "--hclust"), default="complete",
	help="hierarchical clustering method. Choose among the method of the function hclust(). [default=%default]"),

make_option(c("-B", "--base_size"), default=16,
	help="The font base size as defined in ggplot2. [default=%default]"),

make_option(c("-W", "--width"), type="integer", default=7,
	help="Choose the heatmap width in inches. [default=%default]"),

make_option(c("-H", "--height"), type="integer", default=7,
	help="Choose the heatmap height in inches. [default=%default]"),

make_option(c("-o", "--output"), 
	help="Output file name, with the extension. [default=%default]", default="ggdendro.out.pdf"),

make_option(c("-v", "--verbose"), default=FALSE, action="store_true",
	help="Verbose output [default=%default]"),

make_option(c("--debug"), type='integer',
	help="number of lines you want to keep for debugging")

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

# ==========================================
# Function for loading Rdata
# ==========================================

load_obj <- function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}

# ======================
# Plotting variables
# ======================

base_size = opt$base_size
theme_set(theme_grey(base_size))
theme_update(
	axis.ticks=element_blank(),
	axis.ticks.margin = unit(0.01, "inch"),
	axis.ticks.length = unit(0.01, "inch"),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank()
)
theme_update(panel.background = element_blank())



# ===== #
# BEGIN #
# ===== #


# read table
if (opt$verbose) {cat(sprintf("%s: ", Sys.time()), "Reading matrix... ")}
if (opt$input_matrix == "stdin") {
	m = read.table(file("stdin"), h=T)
} else {
	m <- try(load_obj(opt$input_matrix), silent=T)
	if (class(m) == "try-error") {m <- read.table(opt$input_matrix)}
}
if (opt$verbose) {cat("DONE\n")}

set.seed(1)
if (!is.null(opt$debug)) {m = m[sample(nrow(m), opt$debug),]}

# Read palette
if (!is.null(opt$colSide_palette)) {
	colSide_palette_files = strsplit(opt$colSide_palette, ",")[[1]]
	colSide_palette = sapply(colSide_palette_files, function(x)  as.character(read.table(x, h=F, sep="\t", comment.char="%")$V1), simplify=FALSE)
	#if (opt$verbose) {cat("ColSide Palette:", colSide_palette[1], "\n")}
}

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
if (opt$verbose){cat("Removing columns", sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols), "\n")}
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# apply the log10 if needed
if (!is.null(opt$log)) {m <- log(replace(m, is.na(m), 0) + opt$pseudocount, base=eval(parse(text=opt$log)))}


# --------------- Metadata processing -------------

# read metadata
if (opt$verbose) {cat(sprintf("%s: ", Sys.time()), "Read metadata... ")}
if (!is.null(opt$col_metadata)) {
	col_mdata = read.table(opt$col_metadata, h=T, sep="\t", quote="\"", comment.char="")
	# read which fields are needed from the metadata
	if (!is.null(opt$colSide_by)) {colSide_by = strsplit(opt$colSide_by, ",")[[1]]} else {colSide_by = NULL}
	if (!is.null(opt$col_labels) && opt$col_labels != "none") {col_label_fields = strsplit(opt$col_labels,",")[[1]]} else {col_label_fields=NULL}
	
	col_mdata_header = unique(c(opt$merge_col_mdata_on, colSide_by, col_label_fields))
	col_mdata[opt$merge_col_mdata_on] <- gsub(",", ".", col_mdata[,opt$merge_col_mdata_on])
	
	# Select metadata to match the column names (NB: The column Var2 stays)
	df = col_mdata[col_mdata[,opt$merge_col_mdata_on] %in% colnames(m), col_mdata_header]
	colnames(df)[which(colnames(df) == opt$merge_col_mdata_on)] <- "Var2"
}
if (opt$verbose) {cat("DONE\n")}
if (opt$verbose) {print(head(df))}


# ---------------- Dendrogram ----------------------

# COLUMNS

if (opt$vertical) {
	dendro_hjust=0
} else {
	dendro_hjust=1
	dendro_angle=90
}

if (opt$verbose) {cat(sprintf("%s: ", Sys.time()), "Computing distances... ")}
if (opt$dist == "p" || opt$dist =="s") {
	colDist = as.dist(1-cor(m, method=opt$dist, use="p"))
} else {
	colDist = dist(t(m), method=opt$dist)
}
if (opt$verbose) {cat("DONE\n")}

colHC = hclust(colDist, method=opt$hclust)

colHC_data = dendro_data(as.dendrogram(colHC))

col_ggdendro = ggplot(segment(colHC_data))
col_ggdendro = col_ggdendro + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
#col_ggdendro = col_ggdendro + geom_text(data=label(colHC_data), 
#	aes(x=x, y=y, label=label), angle=dendro_angle, hjust=dendro_hjust) 
col_ggdendro = col_ggdendro + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00), "inch")) # top, right, bottom, left
col_ggdendro = col_ggdendro + theme_dendro()
col_ggdendro = col_ggdendro + labs(x=NULL, y=NULL)
col_ggdendro = col_ggdendro + scale_x_continuous(expand=c(0.0, 0.5), labels=NULL) 
col_ggdendro = col_ggdendro + scale_y_continuous(expand=c(0.0, 0.0), labels=NULL)
if (opt$vertical) {
	col_ggdendro = col_ggdendro + coord_flip()
	col_ggdendro = col_ggdendro + scale_y_reverse(expand=c(0.0, 0.0), labels=NULL)
} 
#col_ggdendro



# ---------- Column labels ----------


col_labels_inches = 0
labelData = label(colHC_data)
if (!is.null(opt$col_labels) && opt$col_labels != "none") {
	#labelData$label <- df[,col_label_fields][match(df[,"Var2"], label(colHC_data)[,"label"])]
	labelData$label <- df[,col_label_fields][match(label(colHC_data)[,"label"], df[,"Var2"])]
}


if (is.null(opt$col_labels) || opt$col_labels != "none") {
	gp_labels = ggplot(data=labelData)
	gp_labels = gp_labels + geom_text(aes(x, 0, label=label, hjust=0))
	gp_labels = gp_labels + labs(x=NULL, y=NULL)
	gp_labels = gp_labels + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00), "inch")) # top, right, bottom, left
	gp_labels = gp_labels + scale_x_continuous(expand=c(0.0, 0.5), labels=NULL) 
	gp_labels = gp_labels + scale_y_continuous(expand=c(0.0, 0.0), limits=c(0,1), labels=NULL)
	if (opt$vertical) {
		gp_labels = gp_labels + coord_flip()
	}
	col_labels_inches = max(strwidth(labelData$label, units="in", cex=base_size*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))
}


# -------------------- Column Side Colors ------------

ColSides = list(); ColSide_legends = list()
col_limits = label(colHC_data)[,3]

#if (opt$verbose) {print(col_limits)}

if (!is.null(opt$colSide_by)) {
	if (!is.null(opt$colSide_palette)) {
		if (length(colSide_palette) == 1) {
			colSide_palette = rep(colSide_palette, length(colSide_by))}
		if (length(colSide_palette) >1 && length(colSide_palette) != length(colSide_by)) {
			cat("ERROR: Inconsistent number of column factors and palettes\n"); 
			q(save='no')
		}
    }
	i=1;
	for (colSide in colSide_by) {
		colSide_data = unique(df[c("Var2", colSide)])
		ColSide = ggplot(colSide_data, aes(x=Var2, y="a"))
		ColSide = ColSide + geom_tile(aes_string(fill=colSide))
		ColSide = ColSide + scale_x_discrete(limits = col_limits, labels=NULL, expand=c(0,0))
		ColSide = ColSide + scale_y_discrete(labels=NULL, expand=c(0,0))
		if (!is.null(opt$colSide_palette)) {
			ColSide = ColSide + scale_fill_manual(values=colSide_palette[[i]])
		} else {
			ColSide = ColSide + scale_fill_hue()
		}
#
#		if (!is.null(opt$colSide_palette)) {
#			ColSide = ColSide + scale_fill_manual(values=colSide_palette)
#		}
		ColSide = ColSide + theme(plot.margin=unit(c(0.00, 0.00, 0.00, 0.00),"inch"))
		ColSide = ColSide + labs(x=NULL, y=NULL)
		if (opt$vertical) {
			ColSide = ColSide + coord_flip()
		} 
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

TOT_H = opt$height
TOT_W = opt$width

# >>> Margins <<<

left_margin = 0.1
right_margin = 0.1
vert_inter_space = 0.1
colbar_h = 0.50          # height of each column bar

# >>> Dendro viewport (default) parameters <<<<<

dendro_y = colbar_h*length(ColSides) # + col_labels_inches
dendro_h = TOT_H - 0.05
dendro_w = TOT_W - (right_margin + left_margin)

if (opt$vertical) {
	dendro_y = 0.05
	dendro_h = TOT_H - (right_margin + left_margin)
	dendro_w = TOT_W - 0.05
	dendro_w = dendro_w - col_labels_inches
	dendro_w = dendro_w - (colbar_h+0.02)*length(ColSides)
}

if (!is.null(opt$colSide_by)) {
	legend_width_inch = max(sapply(ColSide_legends, function(x) sum(sapply(x$widths, convertUnit, "in"))))
	legend_height_inch = sum(sapply(ColSide_legends, function(x) sum(sapply(x$heights, convertUnit, "in")))) + 0.1*(length(ColSide_legends)-1)
	dendro_w = dendro_w - (legend_width_inch + vert_inter_space)
}

if (opt$verbose) {cat(sprintf("%s: ", Sys.time()), "Creating viewports... \n")}
if (opt$verbose) {cat('col legend vp.. \n')}

# >>>>> Column side legends viewport <<<<<<<<<<<<<<<<<

if (!is.null(opt$colSide_by)) {
#	legend_width_inch = max(sapply(ColSide_legends, function(x) sum(sapply(x$widths, convertUnit, "in"))))
#	legend_height_inch = sum(sapply(ColSide_legends, function(x) sum(sapply(x$heights, convertUnit, "in")))) + 0.1*(length(ColSide_legends)-1)
	#dendro_w = TOT_W - (legend_width_inch + right_margin + left_margin + vert_inter_space) 
#	dendro_w = dendro_w - (legend_width_inch + vert_inter_space) 
	side_legend_vps = list()
	legend_h = 0
	for (i in 1:length(colSide_by)) {
		side_legend_vp = viewport(
			y = (TOT_H - 0.01) - legend_h,
			x = left_margin + dendro_w,
			h = sum(sapply(ColSide_legends[[i]]$heights, convertUnit, "in")),
			w = legend_width_inch,
			default.units = "inch",
			just = c("left", "top")
		)
		if (opt$vertical) {
			side_legend_vp = viewport(
				y = (TOT_H - 0.01) - legend_h,
				x = left_margin + dendro_w + vert_inter_space + col_labels_inches + (colbar_h+0.02)*length(ColSides),
				h = sum(sapply(ColSide_legends[[i]]$heights, convertUnit, "in")),
				w = legend_width_inch,
				default.units = "inch",
				just = c("left", "top")
			)
		}
		legend_h = sum(sapply(ColSide_legends[[i]]$heights, convertUnit, "in"))
		side_legend_vps[[i]] = side_legend_vp
	}
} else {legend_width_inch = 0}

if (opt$verbose) {cat('col side vp... \n')}

# >>>>>> Column side viewport <<<<<<<<<<<

if (!is.null(opt$colSide_by)){
	ColSide_vps = list(); ColSide_label_vps = list()
	for (i in 1:length(ColSides)) {
		ColSide_vps[[i]] = viewport(
			y = 0.01 + colbar_h*(i-1), 
			x = left_margin,
			h = colbar_h,
			w = dendro_w,
			default.units = "inch",
			just = c("left", "bottom") 
		)
		if (opt$vertical) {
			ColSide_vps[[i]] = viewport(
				y = 0.05,
				x = dendro_w + colbar_h*(i-1),
				h = dendro_h,
				w = colbar_h,
				default.units = "inch",
				just = c("left", "bottom") 
			)
		}
	#	ColSide_label_vps[[i]] = viewport(
	#		y = 0.01 + colbar_h*(i-1),
	#		x = left_margin + dendro_w,
	#		h = colbar_h,
	#		w = as.numeric(strwidth(colSide_by, "inch")),
	#		default.units = "inch",
	#		just = c("left", "bottom")
	#	)
	}
#	if (opt$vertical) {
#		dendro_w = dendro_w - (colbar_h+0.02)*length(ColSides)
#	}
	if (!opt$vertical) {
		dendro_h = TOT_H - (colbar_h+0.02)*length(ColSides)
	}
} 



if (opt$verbose) {cat('dendro vp...\n')}


# >>>>> Column labels viewport <<<<<<<<<<<<<<<<<

if (is.null(opt$col_labels) || opt$col_labels != "none") {
	if (opt$vertical) {
		col_labels_vp = viewport(
			y = 0.05,
			x = dendro_w + vert_inter_space + (colbar_h+0.02)*length(ColSides),
			h = dendro_h,
			w = col_labels_inches,
			default.units = "inch",
			just = c("left", "bottom") 
		)
	}
} 

if (opt$verbose) {cat('col side vp... \n')}


# >>> Column dendrogram viewport <<<<<<<<<

dendro_vp_y = dendro_y # + col_labels_inches
dendro_vp_x = 0.1

colDendro_vp = viewport(
	y = dendro_vp_y, 
	x = left_margin,
	h = dendro_h,
	w = dendro_w,
	default.units = "inch",
	just = c("left", "bottom")
)




# =======================================
# PRINT PLOT
# =======================================

pdf(opt$output, h = TOT_H, w=TOT_W)

#X11(h=total_h, w=total_w)


# Print column side colors
if (!is.null(opt$colSide_by)) {
	for (i in 1:length(ColSide_vps)) {
		print(ColSides[[i]], vp=ColSide_vps[[i]], newpage=FALSE)
#		print(gp_labels, vp=ColSide_vps[[i]], newpage=FALSE)
	#	grid.text(colSide_by[i], x=unit(1,"npc"), just="left", vp=ColSide_label_vps[[i]], gp=gpar(face="bold"))
	}
}

if (opt$verbose) {cat('print col side colors\n')}

## Print column dendrogram
print(col_ggdendro, vp=colDendro_vp, newpage=FALSE)


# Print column side color scales
if (!is.null(opt$colSide_by)) {
	for (i in length(colSide_by):1) {
		all_side_legends = ColSide_legends
		pushViewport(side_legend_vps[[i]]); grid.draw(all_side_legends[[i]]); upViewport()
	}
}

# Print dendrogram labels
if (is.null(opt$col_labels) || opt$col_labels != "none") {
	print(gp_labels, vp=col_labels_vp, newpage=FALSE)
}


dev.off()

file.remove("Rplots.pdf")
q(save="no")
