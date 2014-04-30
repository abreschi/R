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
make_option(c("-x", "--x_axis"), type='integer', help="the index (1-based) of the column you want on the x axis"),
make_option(c("-y", "--y_axis"), type='integer', help="the index (1-based) of the column you want on the y axis"),
make_option(c("-o", "--output_suffix"), help="output filename [default=%default]", default='scatterplot.out.pdf'),
make_option(c("-t", "--type"), help="<tile>, <hex>, <scatter> [default=%default]", default="tile"),
make_option(c("-b", "--binwidth"), help="comma-separated values for binwidth x,y [default=%default]", default="1,1"),
make_option(c("--x_log"), action="store_true", help="x values log10 transformed [default=%default]", default=FALSE),
make_option(c("--y_log"), action="store_true", help="y values log10 transformed [default=%default]", default=FALSE),
make_option(c("--x_psd"), help="pseudocount for x values [default=%default]", default=x_psd, type='double'),
make_option(c("--y_psd"), help="pseudocount for y values [default=%default]", default=y_psd, type='double'),
make_option("--x_title", help="write a title for x axis [default=%default]", default="x_title"),
make_option("--y_title", help="write a title for y axis [default=%default]", default="y_title"),
make_option("--legend_title", help="write a title for the legend [default=%default]", default="count"),

make_option(c("--highlight"), 
	help="a list of element you want to overlay as extra dots"),

make_option(c("--id_col"), default=1,
	help="column with ids"),

make_option("--title", default="", 
	help="write a title for the plot [default=%default]"),

make_option("--diagonal", action="store_true", default=FALSE, 
	help="plot the diagonal [default=%default]"),

make_option(c("-R", "--linear_regression"), action="store_true", default=FALSE, 
	help="plot the regression line [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE, 
	help="verbose output [default=%default]")
)


parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list,
	description = "Plot a density scatterplot"
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


if (opt$x_log) {m[,opt$x_axis] <- m[,opt$x_axis] + opt$x_psd}
if (opt$y_log) {m[,opt$y_axis] <- m[,opt$y_axis] + opt$y_psd}

df = m

# Pearson correlation coefficient
pearson = round(cor(sapply(df[,opt$x_axis], function(x) ifelse(opt$x_log, log10(x), x)), 
	sapply(df[,opt$y_axis], function(x) ifelse(opt$y_log, log10(x), x)), method='p', use='p'), 2)
spearman = round(cor(sapply(df[,opt$x_axis], function(x) ifelse(opt$x_log, log10(x), x)), 
	sapply(df[,opt$y_axis], function(x) ifelse(opt$y_log, log10(x), x)), method='s', use='p'), 2)


# PLOTTING ...

theme_set(theme_bw(base_size=16))

bwidth = as.numeric(strsplit(opt$binwidth, ",")[[1]])
plot_title = sprintf("%s (p_r=%s; s_r=%s)", opt$title, pearson, spearman)
x_col = colnames(df[opt$x_axis])
y_col = colnames(df[opt$y_axis])

# Read the subset of elements you want to highlight
if (!is.null(opt$highlight)) {
	highlight = read.table(opt$highlight, h=F)$V1
	df_h = df[ df[,opt$id_col] %in%  highlight, ]
	if (opt$verbose) {
		print(head(highlight))
		print(head(df_h))
	}
}


countBins <- c(0,1,2,5,10,25,50,75,100,500,Inf)

gp = ggplot(df, aes_string(x=x_col, y=y_col)) 

if (opt$type == 'tile') {
	gp = gp + stat_bin2d(bins=100)
	gp = gp + scale_fill_gradientn(colours=terrain.colors(20), name=opt$legend_title)
	if (!is.null(opt$highlight)) {
		gp = gp + geom_point(data=df_h, aes_string(x=x_col, y=y_col))
	}
}

if (opt$type == 'hex') {
gp = gp + geom_hex(aes(fill=cut(..count.., c(0,1,2,5,10,25,50,75,100,500,Inf))), binwidth=bwidth)
gp = gp + scale_fill_manual('counts', values=terrain.colors(length(countBins))) }

if (opt$type == "scatter") {
gp = gp + geom_point(shape=".")
}


gp = gp + labs(x=opt$x_title, y=opt$y_title, title=plot_title)

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

ggsave(opt$output, h=5, w=6)





q(save='no')

