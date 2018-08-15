#!/usr/bin/env Rscript

# DEFAULT OPTIONS

opt = list()
opt$log10 = FALSE
opt$pseudocount = 1e-04
opt$row_as_variables = FALSE

suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze. Can be stdin"),
make_option(c("-l", "--log10"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-p", "--pseudocount"), type="double", help="specify a pseudocount for the log [default=%default]", default=1e-04),
make_option(c("-m", "--metadata"), help="A list of tsv files with metadata on matrix experiment.\n\t\tThey must be in the format 'file1.tsv,file2.tsv' and contain a key column named 'labExpId'. Can be omitted"),

make_option(c("--merge_mdata_on"), default="labExpId",
	help="[default=%default]"),

#make_option(c("-o", "--output"), help="additional info you want to put in the output file name", default="out"),
make_option(c("-c", "--color_by"), help="choose the fields in the metadata you want to color by", type='character'),

make_option(c("--sort_color"), type='character', 
	help="A field for sorting colors. Can be omitted [default=%default]"),

make_option(c("-s", "--shape_by"), default=NULL, type="character", help="choose the fields in the metadata you want to shape by"),

make_option(c("--no_legend"), action="store_true", default=FALSE,
	help="Do not show the legend [default=%default]"),

make_option(c("-r", "--row_as_variables"), action="store_true", help="select this if you want rows as variables [default=%default]", default=FALSE),
make_option(c("-C", "--princomp"), help="choose the principal components you want to plot. With 3 PC it gives a 3d plot [default='PC1,PC2']", default="PC1,PC2"),

make_option(c("--print_scores"), action="store_true", default=FALSE, 
	help="Output the resuling PCs as a separate file with the extension PCs.tsv [default=%default]"),

make_option(c("--print_loadings"), action="store_true", default=FALSE, 
	help="Output the resulting loadings as a separate file with the extension loadings.tsv [default=%default]"),

make_option(c("--print_lambdas"), action="store_true", default=FALSE,
	help="Output the resulting lambdas (stdev) as a separate file with the extension lambdas.tsv [default=%default]"),

make_option(c("--biplot"), default=FALSE, action="store_true",
	help="If active, the factor of the color is used as grouping factor.
	Centroids are computed and the first <top> loadings are plotted wrt to the two specified components [default=%default]"),

make_option(c("--palette"), default="/users/rg/abreschi/R/palettes/cbbPalette1.15.txt",
	help="File with the color palette [default=%default]"),

make_option(c("--border"), default=FALSE, action="store_true",
	help="Black border to dots [default=%default]"),

make_option(c("--shapes"), 
	help="File with the shapes [default=%default]"),

make_option(c("-L", "--labels"), default=NULL, type="character",
	help="The metadata field with the labels [default=%default]"),

make_option(c("-B", "--base_size"), default=16, type='numeric',
	help="Base font size [default=%default]"),

make_option(c("-H", "--height"), default=7,
	help="Height of the plot in inches [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help="Width of the plot in inches [default=%default]"),

make_option(c("-o", "--output"), default="pca.out",
	help="output file name [default=%default]"),

make_option(c("-v", "--verbose"), action='store_true', default=FALSE,
	help="verbose output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if (opt$verbose) {print(opt)}
##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("grid"))


###############
# BEGIN
##############

# read input tables
inF = opt$input_matrix; if (opt$input_matrix == "stdin") {inF = file("stdin")}
m = read.table(inF, h=T, sep="\t")
if (opt$verbose) {
	cat("Sample of input matrix:\n")
	print(head(m[,1:10]))
}


# Read the color palette
my_palette = read.table(opt$palette, h=F, comment.char="%", sep="\t")$V1

# Read the color ordering
if (is.null(opt$sort_color)) {
	sort_color=NULL
}else{
	sort_color = strsplit(opt$sort_color, ",")[[1]]
}

# Read the shapes
if (!is.null(opt$shapes)) {
	my_shapes = read.table(opt$shapes, h=F, comment.char="%")$V1
}

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

if (opt$verbose) {sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)}

# apply the log if required
if (opt$log10) {m = log10(replace(m, is.na(m), 0) + opt$pseudocount)}

# apply pca
if (opt$row_as_variable) {
m_pca = prcomp(na.omit(m), center=FALSE, scale.=FALSE)} else{
m_pca = prcomp(t(na.omit(m)), center=FALSE, scale.=FALSE)}

# Scale the scores for biplot
#scaledScores = sweep(m_pca$x, 2, m_pca$sdev / sqrt(nrow(m_pca$x)), "/")
scaledScores = m_pca$x

if (opt$verbose) {print(dim(na.omit(m)))}

# HANDLING METADATA

# add metadata to pca results, they should be input in the command line in the future
if (is.null(opt$color_by)) {color_by=NULL
}else{color_by = color_by = strsplit(opt$color_by, ",")[[1]]}
if (is.null(opt$shape_by)) {shape_by=NULL
}else{shape_by = strsplit(opt$shape_by, ",")[[1]]}

# read metadata, one or more table to be merged on labExpId
if (!is.null(opt$metadata)){
	mdata = read.table(opt$metadata, h=T, sep="\t", row.names=NULL, comment.char="", quote="", check.names=F);
	if (opt$merge_mdata_on %in% colnames(mdata)) {
		mdata[,opt$merge_mdata_on] <- gsub("[,-]", ".", mdata[,opt$merge_mdata_on])
	}

	if (opt$verbose) {cat('append metadata...')}
	
	df = merge(as.data.frame(scaledScores),
	unique(mdata[c(color_by, shape_by, opt$merge_mdata_on, opt$labels)]), by.x='row.names', by.y=opt$merge_mdata_on, all.x=T)
	if (opt$verbose) {cat("DONE\n")}
}else{
	df = as.data.frame(scaledScores)
}


if (opt$verbose) {print(head(df))}

#########
# OUTPUT
#########

output_name = opt$output

# Print text outputs if required

# -- principal components --
if (opt$print_scores) {
	write.table(m_pca$x, sprintf("%s.PCs.tsv", output_name), quote=F, sep="\t")
}

# -- loadings --
if (opt$print_loadings) {
	write.table(sweep(m_pca$rotation, 2, m_pca$sdev, "*"), 
		sprintf("%s.loadings.tsv", output_name), quote=F, sep="\t")
}

# -- lambdas --
if (opt$print_lambdas) {
	perc = round(100*m_pca$sdev/sum(m_pca$sdev), 2)
	variances = round(m_pca$sdev^2/sum(m_pca$sdev^2)*100, 2)
	out_df = data.frame(lambda=m_pca$sdev, perc=perc, var_perc=variances)
	write.table(out_df, sprintf("%s.lambdas.tsv", output_name), quote=F, sep="\t")
}

# Read the required components 
prinComp = strsplit(opt$princomp, ",")[[1]]
prinComp_i = as.numeric(gsub("PC", "", prinComp))

# Get a vector with all the variance percentages
variances = round(m_pca$sdev^2/sum(m_pca$sdev^2)*100, 2)

if (opt$biplot) {

	aggrVar = opt$color_by

	# === Centroids ===

	centroids = aggregate (
		df[,which(grepl("^PC", colnames(df)))],
		by=list(df[[aggrVar]]),
		mean
	)
	centroidsM = centroids[,-1]
	
	# === Loadings ===
	
	vecNorm = function(x) {sqrt(sum(x**2))}
	
	scaledLoadings = sweep(m_pca$rotation, 2, m_pca$sdev, "*")
	
	#centroidsNorm = apply(centroidsM[,prinComp], 1, vecNorm)         # DIM: number of levels x 1
	#loadingsNorm = apply(scaledLoadings[,prinComp], 1, vecNorm)      # DIM: number of variables x 1
	centroidsNorm = apply(centroidsM, 1, vecNorm)         # DIM: number of levels x 1
	
	#cosine = ( scaledLoadings[,prinComp] %*% t(centroidsM[,prinComp]) ) / (loadingsNorm %*% t(centroidsNorm))
	cosine = scaledLoadings %*% t(centroidsM/centroidsNorm) 
	cosine = setNames(data.frame(cosine),  centroids[,1])

	closest = setNames(melt(apply(1-cosine, 2, rank)), c("variable", aggrVar, "rank"))
	write.table( cosine, file=sprintf("%s.cosine.tsv", opt$output), quote=F, sep="\t");

	closest_df = data.frame(merge(closest, scaledLoadings, 
		by.x="variable", by.y="row.names"), check.names=F)

	
}



#############
# PLOT
#############

# plot parameters
pts = 5

l_col = opt$labels
base_size = opt$base_size

theme_set(theme_bw(base_size = base_size))
theme_update(legend.text=element_text(size=0.9*base_size),
	legend.key.size=unit(0.9*base_size, "points"),
	legend.key = element_blank()
)

top = 30


# Open device for plotting
pdf(sprintf("%s.pdf", output_name), w=opt$width, h=opt$height)

if (length(prinComp) == 2){

	geom_params = list()
	geom_params$size = pts
#	geom_params$alpha = opt$alpha
	
	mapping = list()
	mapping <- modifyList(mapping, aes_string(x=prinComp[1], y=prinComp[2]))

	if (!is.null(opt$color_by)) {
		gp_color_by=interaction(df[color_by])
		if (!is.null(opt$sort_color)) {
			gp_color_by = factor(gp_color_by, levels=sort_color)
		}
		mapping = modifyList(mapping, aes_string(color=gp_color_by, order=gp_color_by))
	} else {
		gp_color_by=NULL
	}
	
	if (!is.null(opt$shape_by)) {
		gp_shape_by=interaction(df[shape_by])
		if (!is.null(opt$sort_shape)) {
			gp_shape_by = factor(gp_shape_by, levels=sort_shape)
		}
		mapping = modifyList(mapping, aes_string(shape=gp_shape_by, order=gp_shape_by))
	} else {
		gp_shape_by=NULL
	}

#	if (!is.na(opt$shape_by)) {gp_shape_by=interaction(df[shape_by]);
#	gp_shape_by <- factor(gp_shape_by, levels=sort(levels(gp_shape_by)))
#	mapping = modifyList(mapping, aes_string(shape=S_col))
	
	class(mapping) <- "uneval"
	
	pointLayer <- layer(
		geom = "point",
	#	geom_params = geom_params,
		params = geom_params,
		mapping = mapping,
		stat = "identity",
		position = "identity"
	)
	



	# plotting...
	gp = ggplot(df, aes_string(x=prinComp[1],y=prinComp[2]));

	if (opt$biplot) {
		gp = gp + geom_point(data=centroids, aes_string(
			x=prinComp[1], y=prinComp[2], color="Group.1"), shape=8, size=7)
		gp = gp + geom_segment( 
			data=subset(closest_df, rank <= top), 
			aes_string(x=0, y=0, xend=prinComp[1], yend=prinComp[2], 
				color=sprintf("`%s`", opt$color_by))
		)
	}


	if (opt$border) {
		if (!is.null(opt$shape_by)) {
		gp = gp + geom_point(aes(shape=gp_shape_by), col='black', size=pts+1.0);
		} else {
		gp =  gp + geom_point(col="black", size=pts+1.0)
		}
	}

	gp = gp + pointLayer

#	gp = gp + geom_point(aes(color=gp_color_by))
#	gp = gp + geom_point(aes(col=gp_color_by, shape=gp_shape_by), size=pts);
#
	gp = gp + labs(title="");
	gp = gp + labs(x=sprintf('%s (%s%%)', prinComp[1], variances[prinComp_i[1]]));
	gp = gp + labs(y=sprintf('%s (%s%%)', prinComp[2], variances[prinComp_i[2]]));

	gp = gp + scale_color_manual(name=opt$color_by, values=my_palette)
	if (!is.null(opt$shapes)) {
		gp = gp + scale_shape_manual(name=opt$shape_by, values=my_shapes);
	}
	if (opt$no_legend) {
		gp = gp + guides(shape=FALSE, color=FALSE)
	
}
	if (!is.null(opt$labels)) {
		gp = gp + geom_text(aes_string(label=l_col), size=pts)
	}

	gp
} 




# --------------------
#
# 3d scatterplot
#
# --------------------


if (length(prinComp) == 3) {

suppressPackageStartupMessages(library(scatterplot3d))

par(xpd=NA, omi=c(0.5, 0.5, 0.5, 1.0))

if (!is.na(opt$color_by)) {gp_color=my_palette[interaction(df[color_by])]} else {gp_color="black"}
if (!is.null(opt$shape_by)) {gp_shape_by=interaction(df[shape_by]);
gp_shape_by <- factor(gp_shape_by, levels=sort(intersect(levels(gp_shape_by), gp_shape_by))); gp_shape=my_shapes[gp_shape_by]} else {gp_shape_by=NULL}

plot3d = scatterplot3d(df[prinComp], 
	color = gp_color,
	pch = gp_shape,
	xlab = sprintf('%s (%s%%)', prinComp[1], variances[prinComp_i[1]]),
	ylab = sprintf('%s (%s%%)', prinComp[2], variances[prinComp_i[2]]),
	zlab = sprintf('%s (%s%%)', prinComp[3], variances[prinComp_i[3]]),
	cex.symbols = 1.5,
	lab = c(5,4)
)

# !!! To be removed after the mouse paper !!!
#i=0; for(sample in interaction(df[color_by])) {
#i=i+1; plot3d$points3d(subset(df, General_category == sample, select=prinComp), type='l', col=gp_color[i])}

if (!is.na(opt$color_by)) {
	legend(
		x = log(max(df[prinComp[1]])) + 3,
#		x = 5,
		y = 5.5,
		legend = levels(interaction(df[color_by])), 
		fill = my_palette[seq_along(levels(interaction(df[color_by])))]
	)
}

if (!is.na(opt$shape_by)) {
	legend(
#		x = -log(abs(min(df[prinComp[1]]))) - 1.5, 
		x = -3,
		y = 6, 
#		y = 7.2,
		legend = levels(gp_shape_by), 
		pch = my_shapes[seq_along(levels(gp_shape_by))]
		)
#	legend(-log(abs(min(df[prinComp[1]])))+1.5,7.2,levels(gp_shape_by), 
#	pch=shapes[seq_along(levels(gp_shape_by))])
}
}


dev.off()
q(save='no')

