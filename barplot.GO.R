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

make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

make_option(c("-f", "--fill"), default="dodgerblue",
	help="choose the color which you want to fill the histogram with [default=%default]"),

make_option(c("-F", "--fill_by"), type='integer', default=NULL, 
	help="the column with the factor for the color [default=%default]"),

make_option(c("-P", "--palette"), default=NULL, 
	help="palette file [default=%default]"),

make_option(c("-H", "--height"), default=7,
	help="Height of the plot in inches [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help="Width of the plot in inches [default=%default]"),

make_option(c("-B", "--base_size"), default=20,
	help="Base size [default=%default]"),

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

# Read palette file
if (!is.null(opt$palette)) {
	palette = read.table(opt$palette, h=F, sep="\t", comment.char="")$V1
}


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

theme_set(theme_bw(base_size=opt$base_size))
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

# specify fill column
if (!is.null(opt$fill_by)) {
	F_col = colnames(df)[opt$fill_by]
	mapping <- aes_string(fill=F_col)
} else {
	mapping = NULL
}

# define histogram layer 
histLayer <- layer(
    geom = "bar",
#    geom_params = geom_params,
	params = geom_params,
	position = "identity",
	mapping = mapping,
    stat = "identity"
)


# start the plot
gp = ggplot(df, aes_string(x = x_col, y = y_col))
gp = gp + histLayer
gp = gp + coord_flip()
gp = gp + labs(y="-log10(pvalue)", x=NULL)
if (!is.null(opt$palette)) {
	gp = gp + scale_fill_manual(values=palette)
}
ggsave(opt$output, h=opt$height, w=opt$width, title=opt$output)

q(save='no')
