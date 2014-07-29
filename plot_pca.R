#!/usr/bin/env Rscript

# DEFAULT OPTIONS

opt = list()
opt$log10 = FALSE
opt$pseudocount = 1e-04
opt$row_as_variables = FALSE

##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library("reshape"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)
#custom_palette <- rgb(matrix(c(0,0,0,0,73,73,0,146,146,255,109,182,255,182,119,73,0,146,0,109,219,182,109,
#255,109,182,255,182,219,255,146,0,0,146,73,0,219,209,0,36,255,36,255,255,109), ncol=3, byrow=T), max=255)
#cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#my_palette <- custom_palette

##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze. Can be stdin"),
make_option(c("-l", "--log10"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-p", "--pseudocount"), type="double", help="specify a pseudocount for the log [default=%default]", default=1e-04),
make_option(c("-m", "--metadata"), help="A list of tsv files with metadata on matrix experiment.\n\t\tThey must be in the format 'file1.tsv,file2.tsv' and contain a key column named 'labExpId'. Can be omitted"),
#make_option(c("-o", "--output"), help="additional info you want to put in the output file name", default="out"),
make_option(c("-c", "--color_by"), help="choose the fields in the metadata you want to color by", type='character'),

make_option(c("-s", "--shape_by"), help="choose the fields in the metadata you want to shape by", type='character', default=NA),

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

make_option(c("--palette"), default="/users/rg/abreschi/R/palettes/cbbPalette1.15.txt",
	help="File with the color palette [default=%default]"),

make_option(c("-H", "--height"), default=7,
	help="Height of the plot in inches [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help="Width of the plot in inches [default=%default]"),

make_option(c("-o", "--output"), default="pca.out",
	help="output file name"),

make_option(c("-v", "--verbose"), action='store_true', default=FALSE,
	help="verbose output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if (opt$verbose) {print(opt)}


###############
# BEGIN
##############

# read input tables
if (opt$input_matrix == "stdin") {
	m = read.table(file("stdin"), h=T)
} else {
	m = read.table(opt$input_matrix, h=T)
}

# Read the color palette
my_palette = read.table(opt$palette, h=F, comment.char="%")$V1

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


if (opt$verbose) {print(dim(na.omit(m)))}

# HANDLING METADATA

# add metadata to pca results, they should be input in the command line in the future
if (is.null(opt$color_by)) {color_by=NULL
}else{color_by = color_by = strsplit(opt$color_by, ",")[[1]]}
if (is.na(opt$shape_by)) {shape_by=NULL
}else{shape_by = strsplit(opt$shape_by, ",")[[1]]}

# read metadata, one or more table to be merged on labExpId
if (!is.null(opt$metadata)){
	metadata = strsplit(opt$metadata, ",")[[1]]
	for (i in seq_along(metadata)) {
		mdata = read.table(metadata[i], h=T, sep="\t", row.names=NULL);
		if ('labExpId' %in% colnames(mdata)) {
			mdata$labExpId <- gsub(",", ".", mdata[,"labExpId"])
		}
		if ( i==1 ) {
			new_mdata = mdata
		}else{
		new_mdata = merge(mdata, new_mdata, by=c('cell'))
		}
	}

	cat('append metadata...')
	
	df = merge(as.data.frame(m_pca$x),
	unique(new_mdata[c(color_by, shape_by, 'labExpId')]), by.x='row.names', by.y='labExpId', all.x=T)
}else{
	df = as.data.frame(m_pca$x)
}




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
	write.table(m_pca$rotation, sprintf("%s.loadings.tsv", output_name), quote=F, sep="\t")
}

# -- lambdas --
if (opt$print_lambdas) {
	write.table(m_pca$sdev, sprintf("%s.lambdas.tsv", output_name), quote=F, sep="\t")
}

# Read the required components 
prinComp = strsplit(opt$princomp, ",")[[1]]
prinComp_i = as.numeric(gsub("PC", "", prinComp))

# Get a vector with all the variance percentages
variances = round(m_pca$sdev^2/sum(m_pca$sdev^2)*100, 2)


#############
# PLOT
#############

# plot parameters
pts = 5
#shapes = c(16, 15, 0, 14:1, 18, 17)
shapes = c(15, 0)
theme_set(theme_bw(base_size=16))
theme_update(legend.key = element_blank())


# Open device for plotting
pdf(sprintf("%s.pdf", output_name), w=opt$width, h=opt$height)

if (length(prinComp) == 2){
	# plotting...
	gp = ggplot(df, aes_string(x=prinComp[1],y=prinComp[2]));
	if (!is.null(opt$color_by)) {gp_color_by=interaction(df[color_by])} else {gp_color_by=NULL}
	if (!is.na(opt$shape_by)) {gp_shape_by=interaction(df[shape_by]);
	gp_shape_by <- factor(gp_shape_by, levels=sort(levels(gp_shape_by)))} else {gp_shape_by=NULL}
	gp = gp + geom_point(aes(col=gp_color_by, shape=gp_shape_by), size=pts);
	gp = gp + labs(title="");
	gp = gp + labs(x=sprintf('%s (%s%%)', prinComp[1], variances[prinComp_i[1]]));
	gp = gp + labs(y=sprintf('%s (%s%%)', prinComp[2], variances[prinComp_i[2]]));
	gp = gp + scale_color_manual(name=opt$color_by, values=my_palette)
	gp = gp + scale_shape_manual(name='Sample', values=shapes);
	if (opt$no_legend) {
		gp = gp + guides(shape=FALSE, color=FALSE)
	}
	gp
} 




# --------------------
#
# 3d scatterplot
#
# --------------------

#print(head(df))

#if (length(prinComp) == 3) {
#
#	suppressPackageStartupMessages(library(ggtern))
#
#	gp = ggtern(df, aes_string(x=prinComp[1], y=prinComp[2], z=prinComp[3]));
#	gp = gp + geom_point()
##	if (!is.na(opt$color_by)) {gp_color_by=interaction(df[color_by])} else {gp_color_by=NULL}
##	if (!is.na(opt$shape_by)) {gp_shape_by=interaction(df[shape_by]);
##	gp_shape_by <- factor(gp_shape_by, levels=sort(levels(gp_shape_by)))} else {gp_shape_by=NULL}
##	gp = gp + geom_point(aes(col=gp_color_by, shape=gp_shape_by), size=pts);
##	gp = gp + labs(title="");
##	gp = gp + labs(x=sprintf('%s (%s%%)', prinComp[1], variances[prinComp_i[1]]));
##	gp = gp + labs(y=sprintf('%s (%s%%)', prinComp[2], variances[prinComp_i[2]]));
##	gp = gp + scale_color_manual(name=opt$color_by, values=my_palette)
##	gp = gp + scale_shape_manual(name='Sample', values=shapes);
#	gp
##}

if (length(prinComp) == 3) {

suppressPackageStartupMessages(library(scatterplot3d))

par(xpd=NA, omi=c(0.5, 0.5, 0.5, 1.0))

if (!is.na(opt$color_by)) {gp_color=my_palette[interaction(df[color_by])]} else {gp_color="black"}
if (!is.na(opt$shape_by)) {gp_shape_by=interaction(df[shape_by]);
gp_shape_by <- factor(gp_shape_by, levels=sort(intersect(levels(gp_shape_by), gp_shape_by))); gp_shape=shapes[gp_shape_by]} else {gp_shape_by=NULL}

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
		pch = shapes[seq_along(levels(gp_shape_by))]
		)
#	legend(-log(abs(min(df[prinComp[1]])))+1.5,7.2,levels(gp_shape_by), 
#	pch=shapes[seq_along(levels(gp_shape_by))])
}
}


dev.off()
q(save='no')

