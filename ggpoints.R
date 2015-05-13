#!/usr/bin/env Rscript

options(stringsAsFactors=F)
x_psd = 1e-03
y_psd = 1e-03




##################
# OPTION PARSING
##################
suppressPackageStartupMessages(library("optparse"))


option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin",
	help="the matrix you want to analyze. \"stdin\" for stdin [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE, help="The file has header [default=%default]"),
make_option(c("--xy"), type='character', default="1,2",
	help="the indeces (1-based) of the columns you want on the x axis and on the y axis, comma-separated [default=%default]"),

make_option(c("-C", "--color_by"), type="integer", 
	help="Index of the column by which to color the dots [default=%default]"),

make_option(c("-P", "--palette"), type="character",
	help="File with the palette for colors. Leave empty for color_hue"),

make_option(c("-o", "--output_suffix"), help="output filename [default=%default]", default='scatterplot.out.pdf'),
#make_option(c("-b", "--binwidth"), help="comma-separated values for binwidth x,y [default=%default]", default="1,1"),
make_option(c("--x_log"), action="store_true", help="x values log10 transformed [default=%default]", default=FALSE),
make_option(c("--y_log"), action="store_true", help="y values log10 transformed [default=%default]", default=FALSE),
make_option(c("--x_psd"), help="pseudocount for x values [default=%default]", default=x_psd, type='double'),
make_option(c("--y_psd"), help="pseudocount for y values [default=%default]", default=y_psd, type='double'),
make_option("--x_title", help="write a title for x axis [default=%default]", default="x_title"),
make_option("--y_title", help="write a title for y axis [default=%default]", default="y_title"),
make_option("--legend_title", help="write a title for the legend [default=%default]", default="count"),

make_option("--title", default="", 
	help="write a title for the plot [default=%default]"),

make_option("--diagonal", action="store_true", default=FALSE, 
	help="plot the diagonal [default=%default]"),

make_option(c("-R", "--linear_regression"), action="store_true", default=FALSE, 
	help="plot the regression line [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help="Width for the plot in inches [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE, 
	help="verbose output [default=%default]")
)


parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list,
	description = "Plot a scatterplot"
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

##------------
## LIBRARIES
##------------ 

if (opt$verbose) {cat("Loading libraries... ")}
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}




###################
# BEGIN           #
###################

if (opt$input_matrix == "stdin") {
	m = read.table(file("stdin"), h=opt$header)
} else {
	m = read.table(opt$input_matrix, h=opt$header)
}

# Read color palette
if (!is.null(opt$palette)) {palette = as.character(read.table(opt$palette, h=F, comment.char="%")$V1)}

# Read the columns indeces
axes = strsplit(opt$xy, ",")[[1]]
x_axis = as.integer(axes[1])
y_axis = as.integer(axes[2])

if (opt$x_log) {m[,x_axis] <- m[,x_axis] + opt$x_psd}
if (opt$y_log) {m[,y_axis] <- m[,y_axis] + opt$y_psd}

df = m


# Pearson correlation coefficient
pearson = round(cor(sapply(df[,x_axis], function(x) ifelse(opt$x_log, log10(x), x)), 
	sapply(df[,y_axis], function(x) ifelse(opt$y_log, log10(x), x)), method='p', use='p'), 2)
spearman = round(cor(sapply(df[,x_axis], function(x) ifelse(opt$x_log, log10(x), x)), 
	sapply(df[,y_axis], function(x) ifelse(opt$y_log, log10(x), x)), method='s', use='p'), 2)


# PLOTTING ...

theme_set(theme_bw(base_size=16))
theme_update(
	legend.key=element_blank(),
	title = element_text(vjust=1),
	axis.title.x = element_text(vjust=0)
)

plot_title = sprintf("%s (p_r=%s; s_r=%s)", opt$title, pearson, spearman)
x_col = colnames(df[x_axis])
y_col = colnames(df[y_axis])
C_col = colnames(df)[opt$color_by]

gp = ggplot(df, aes_string(x=x_col, y=y_col)) 
if (!is.null(opt$color_by)) {
	gp = gp + geom_point(aes_string(color=C_col), size=1)
} else {
	gp = gp + geom_point(, size=1)
}


if (!is.null(opt$palette)) {
	if (is.character(df[,C_col]) | is.factor(df[,C_col])) {
		gp = gp + scale_color_manual(values=palette)
	} else {
		gp = gp + scale_color_gradientn(colours=palette)
	}
}
gp = gp + labs(x=opt$x_title, y=opt$y_title, title=plot_title)
gp = gp + guides(colour = guide_legend(override.aes = list(size=5))) 

# Add the diagonal line
if (opt$diagonal) {
	gp = gp + geom_abline(intercept=0, slope=1, color="grey")
}

# Add the regression line

if (opt$linear_regression) {
	
	if (opt$verbose) {print(head(df))}
	if (opt$x_log) {
		x_col = sprintf("log10(%s)", x_col)
	}
	if (opt$y_log) {
		y_col = sprintf("log10(%s)", y_col)
	}
	formula = as.formula(sprintf("%s~%s", y_col, x_col))
	coeff = lm(formula, df)$coefficients
	gp = gp + geom_abline(intercept=coeff[1], slope=coeff[2])
}


# Change to log scale

if (opt$x_log) {gp = gp + scale_x_log10()}
if (opt$y_log) {gp = gp + scale_y_log10()}

ggsave(opt$output, h=5, w=opt$width)
	




q(save='no')
