#!/usr/bin/env Rscript


options(stringsAsFactors=FALSE)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="File or stdin [default=%default]"),

make_option(c("-o", "--output"), default="gghist.out.pdf",
	help="Output file name [default=%default]"),

make_option(c("-x", "--x_axis"), default=1,
	help="Index of the column with values, or labels if you already have counts [default=%default]"),

make_option(c("-y", "--y_axis"), default=NULL, type="integer",
	help="Index of the column with values, in case x provides counts. This will plot identity. Leave empty for default histogram [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE,
	help="Use this if the input has a header [default=%default]"),

make_option(c("--position"), default='dodge',
	help="Position for histogram [default=%default]"),

make_option(c("--scale_x_log10"), action="store_true", default=FALSE,
	help="log10-transform x scale [default=%default]"),

make_option(c("--scale_y_log10"), action="store_true", default=FALSE,
	help="log10-transform y scale [default=%default]"),

make_option(c("--y_title"), type="character", default="count",
	help="Title for the y axis [default=%default]"),

make_option(c("--x_title"), type="character", default=NULL,
	help="Title for the x axis [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]"),

make_option(c("-f", "--fill"), default="aquamarine",
	help="choose the color which you want to fill the histogram with"),

make_option(c("-c", "--color"), default="white",
	help="choose the color which you want to contour the histogram with"),

make_option(c("-F", "--fill_by"), type='numeric',
	help="the column index with the factor to fill by. Leave empty for no factor."),

make_option(c("-P", "--palette"), 
        help='File with colors for the lines. Leave empty to use even color spacing'),

make_option(c("--facet_by"), type='numeric',
	help="the column index with the factor to facet by. Leave empty for no factor."),

make_option(c("-W", "--width"), default=7,
	help="width of the plot in inches. [default=%default]"),

make_option(c("-b", "--binwidth"), type="double", 
	help="Specify binwidth. Leave empty for default")

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
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


# ======== #
# BEGIN    #
# ======== #


# Read data
if (opt$input == "stdin") {input=file("stdin")} else {input=opt$input}
m = read.table(input, sep="\t", h=opt$header) 

df = m

# Read facet
if (!is.null(opt$facet_by)) {facet_formula = as.formula(sprintf("~%s", colnames(df)[opt$facet]))}

# Read columns
x_col = colnames(df)[opt$x_axis]
if (!is.null(opt$y_axis)) {y_col = colnames(df)[opt$y_axis]}
if (!is.null(opt$fill_by)) {F_col = colnames(df)[opt$fill_by]}

# Read palette
if (!is.null(opt$palette)) {
	palette = as.character(read.table(opt$palette, h=F, comment.char="%")$V1)
}


df[,x_col] <- gsub("\\\\n", "\n", df[,x_col])

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

# Stat parameters 
stat = ifelse(is.null(opt$y_axis), "bin", "identity")

stat_params = list(
	right=TRUE, 
	include.lowest=TRUE
)

# specify binwidth
if (!is.null(opt$binwidth)) {
	stat_params$binwidth = opt$binwidth
}

geom_params = list()

if (is.null(opt$fill_by)) {
	geom_params$fill = opt$fill
	geom_params$color = opt$color
}

# specify fill column
if (!is.null(opt$fill_by)) {
	mapping <- aes_string(fill=F_col)
} else {
	mapping = NULL
}

# define histogram layer 
histLayer <- layer(
    geom = "bar",
    geom_params = geom_params,
	position = opt$position,
	mapping = mapping,
    stat = stat,
    stat_params = stat_params
)


# start the plot
if (is.null(opt$y_axis)) {
	gp = ggplot(df, aes_string(x = x_col))
} else {
	gp = ggplot(df, aes_string(x = x_col, y = y_col))
}

gp = gp + histLayer

if (!is.character(df[,x_col])) {
	avg = mean(df[,x_col], na.rm=TRUE)
	med = median(df[,x_col], na.rm=TRUE)
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

if (!is.null(opt$x_title)) {gp = gp + labs(x=opt$x_title)}

gp = gp + labs(y=opt$y_title)

ggsave(opt$output, h=5, w=opt$width, title=opt$output)

# EXIT
quit(save='no')
