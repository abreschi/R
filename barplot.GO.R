#!/usr/bin/env Rscript


options(stringsAsFactors=FALSE)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin [default=%default]"),

make_option(c("-o", "--output"), default="barplot.GO.pdf",
	help="Output file name [default=%default]"),

#make_option(c("-x", "--x_axis"), default=1,
#	help="Index of the column with values, or labels if you already have counts [default=%default]"),
#
#make_option(c("-y", "--y_axis"), default=NULL, type="integer",
#	help="Index of the column with values, in case x provides counts. This will plot identity. Leave empty for default histogram [default=%default]"),
#
make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

#make_option(c("--position"), default='dodge',
#	help="Position for histogram [default=%default]"),
#
#make_option(c("--scale_x_log10"), action="store_true", default=FALSE,
#	help="log10-transform x scale [default=%default]"),
#
#make_option(c("--scale_y_log10"), action="store_true", default=FALSE,
#	help="log10-transform y scale [default=%default]"),
#
#make_option(c("--y_title"), type="character", default="count",
#	help="Title for the y axis [default=%default]"),
#
#make_option(c("--x_title"), type="character", default="",
#	help="Title for the x axis [default=%default]"),
#
make_option(c("-f", "--fill"), default="dodgerbluee",
	help="choose the color which you want to fill the histogram with"),

#make_option(c("-c", "--color"), default="white",
#	help="choose the color which you want to contour the histogram with"),
#
make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")

)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list,
	description = "Reads the values on the first column and outputs a histogram"
)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

#------------
# LIBRARIES
#------------ 

if (opt$verbose) {cat("Loading libraries... ")}
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


# ======== #
# BEGIN    #
# ======== #


# Read data
if (opt$input == "stdin") {input=file("stdin")} else {input=opt$input}
m = read.table(input, h=opt$header, sep="\t") 

df = m

# Read columns
y_col = 2
y_col = colnames(df)[y_col]
df[y_col] <- -log10(df[y_col])
opt$x_axis <- 7
x_col = colnames(df)[opt$x_axis]

# Insert newlines in labels if they are too long
max_nchar = 25
toNL = which(sapply(as.character(df[,x_col]), nchar) > max_nchar)
nbSpaces = str_count(as.character(df[,x_col]), " ") 
toNL = intersect(which(nbSpaces > 1), toNL)
#print(df)
matches = gregexpr(" ", as.character(df[, x_col]), fixed=TRUE)
pos = as.numeric(sapply(matches, function(x) {i=floor((length(as.numeric(x))+1)/2); x[i]}))
subst = sapply(1:length(df[,x_col]), function(i) {x=as.character(df[i,x_col]); substring(x, pos[i])<- "\n"; x }) 
df[toNL, x_col] <- subst[toNL]

# Sort the labels by p-value
levs = df[order(df[,y_col]),x_col]
df[,x_col] <- factor(df[,x_col], levels=levs)

#================
# GGPLOT
#================

theme_set(theme_bw(base_size=20))
theme_update(
	axis.text.x=element_text(angle=45, hjust=1, vjust=1),
	legend.key = element_rect(color='white'),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank()
)

#opt$fill = "dodgerblue"
opt$color = "black"

geom_params = list()

if (is.null(opt$fill_by)) {
	geom_params$fill = opt$fill
	geom_params$color = opt$color
}

## specify fill column
#if (!is.null(opt$fill_by)) {
#	mapping <- aes_string(fill=F_col)
#} else {
#	mapping = NULL
#}
mapping = NULL

# define histogram layer 
histLayer <- layer(
    geom = "bar",
    geom_params = geom_params,
	position = "identity",
	mapping = mapping,
    stat = "identity"
)

opt$width = 7

# start the plot
gp = ggplot(df, aes_string(x = x_col, y = y_col))
gp = gp + histLayer
gp = gp + coord_flip()
gp = gp + labs(y="-log10(pv)", x=NULL)
ggsave(opt$output, h=5, w=opt$width, title=opt$output)

q(save='no')

if (!is.character(df[,1])) {
	avg = mean(df[,1], na.rm=TRUE)
	med = median(df[,1], na.rm=TRUE)
	gp = gp + geom_point(aes(x=avg, y=0), size=2)
	gp = gp + geom_vline(xintercept=med, linetype=2)
}

# Color scale
if (!is.null(opt$fill_by)) {
	if (!is.null(opt$palette)) {
		gp = gp + scale_fill_manual(values=palette)
	} else {
		gp = gp + scale_fill_hue()
	}
}

if (!is.null(opt$facet_by)) {
	gp = gp + facet_wrap(facet_formula)
}

if (opt$scale_x_log10) {gp = gp + scale_x_log10()}
if (opt$scale_y_log10) {gp = gp + scale_y_log10()}

gp = gp + labs(y=opt$y_title, x=opt$x_title)

ggsave(opt$output, h=5, w=opt$width, title=opt$output)

# EXIT
quit(save='no')
