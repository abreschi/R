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

make_option(c("-S", "--shape_by"), type="integer", 
	help="Index of the column by which to shape the dots [default=%default]"),

make_option(c("-f", "--facet_by"), type="integer", default=NULL,
	help="Index of the column by faceting [default=%default]"),

make_option(c("--facet_scale"), default="fixed",
	help="Scale of the facet: < fixed | free_y | free_x | free > [default=%default]"),

make_option(c("--facet_nrow"), type="integer", default=NULL,
	help="Number of rows for faceting [default=%default]"),

make_option(c("-a", "--alpha"), default=1,
	help="Transparency value [default=%default]"),

make_option(c("-s", "--size"), default=1,
	help="Size of the points [default=%default]"),

make_option(c("-o", "--output_suffix"), help="output filename [default=%default]", default='scatterplot.out.pdf'),
#make_option(c("-b", "--binwidth"), help="comma-separated values for binwidth x,y [default=%default]", default="1,1"),
make_option(c("--x_log"), action="store_true", help="x values log10 transformed [default=%default]", default=FALSE),
make_option(c("--y_log"), action="store_true", help="y values log10 transformed [default=%default]", default=FALSE),
make_option(c("--x_psd"), help="pseudocount for x values [default=%default]", default=x_psd, type='double'),
make_option(c("--y_psd"), help="pseudocount for y values [default=%default]", default=y_psd, type='double'),
make_option("--x_title", help="write a title for x axis [default=%default]", default="x_title"),
make_option("--y_title", help="write a title for y axis [default=%default]", default="y_title"),

make_option(c("--x_limits"), 
	help="Specify limits for the x-axis scale, e.g. \"\\-1,1\". Escape character for negative numbers [default=%default]"),

make_option(c("--y_limits"), 
	help="Specify limits for the y-ayis scale, e.g. \"\\-1,1\". Escape character for negative numbers [default=%default]"),

make_option("--legend_title", help="write a title for the legend [default=%default]", default="count"),

make_option("--title", default="", 
	help="write a title for the plot [default=%default]"),

make_option("--abline", default=NULL, type="character",
	help="A list of pairs of intercepts and slopes, e.g. i1,s1;i2,s2 [default=%default]"),

make_option("--diagonal", action="store_true", default=FALSE, 
	help="plot the diagonal [default=%default]"),

make_option("--off_diagonal", action="store_true", default=FALSE,
	help="plot the off-diagonal [default=%default]"),

make_option(c("-R", "--linear_regression"), action="store_true", default=FALSE, 
	help="plot the regression line [default=%default]"),

make_option(c("-W", "--width"), default=7,
	help="Width for the plot in inches [default=%default]"),

make_option(c("-H", "--height"), default=7,
	help="Height for the plot in inches [default=%default]"),

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

if(opt$input_matrix == "stdin") {inF = file("stdin")} else {inF=opt$input_matrix}
m = read.table(inF, h=opt$header, sep="\t")

# Read color palette
if (!is.null(opt$palette)) {palette = as.character(read.table(opt$palette, h=F, comment.char="%", sep="\t")$V1)}

# Read the columns indeces
axes = strsplit(opt$xy, ",")[[1]]
x_axis = as.integer(axes[1])
y_axis = as.integer(axes[2])

if (opt$x_log) {m[,x_axis] <- m[,x_axis] + opt$x_psd}
if (opt$y_log) {m[,y_axis] <- m[,y_axis] + opt$y_psd}

df = m

plot_title = opt$title

# Pearson correlation coefficient
if (all(sapply(df[, c(x_axis, y_axis)], class) == "numeric")) {
	pearson = round(cor(sapply(df[,x_axis], function(x) ifelse(opt$x_log, log10(x), x)), 
		sapply(df[,y_axis], function(x) ifelse(opt$y_log, log10(x), x)), method='p', use='p'), 2)
	spearman = round(cor(sapply(df[,x_axis], function(x) ifelse(opt$x_log, log10(x), x)), 
		sapply(df[,y_axis], function(x) ifelse(opt$y_log, log10(x), x)), method='s', use='p'), 2)
	plot_title = sprintf("%s (p_r=%s; s_r=%s)", opt$title, pearson, spearman)
}

# PLOTTING ...

theme_set(theme_bw(base_size=16))
theme_update(
	legend.key=element_blank(),
	title = element_text(vjust=1),
	axis.title.x = element_text(vjust=0)
)


x_col = colnames(df[x_axis])
y_col = colnames(df[y_axis])
C_col = colnames(df)[opt$color_by]
S_col = colnames(df)[opt$shape_by]


if (!is.null(opt$color_by)) {
	df = df[order(df[,C_col]),]
}


geom_params = list()
geom_params$size = opt$size
geom_params$alpha = opt$alpha


mapping = list()
mapping <- modifyList(mapping, aes_string(x=x_col, y=y_col))
if (!is.null(opt$color_by)) {
	mapping = modifyList(mapping, aes_string(color=C_col))
}

#print(head(df[,C_col], 20))
#print(head(order(order(df[,C_col])), 20))

if (!is.null(opt$shape_by)) {
	mapping = modifyList(mapping, aes_string(shape=S_col))
}
class(mapping) <- "uneval"

pointLayer <- layer(
	geom = "point",
#	geom_params = geom_params,
	params = geom_params,
	mapping = mapping,
	stat = "identity",
	position = "identity"
)

gp = ggplot(df) + pointLayer


# Palette

if (!is.null(opt$palette)) {
	if (is.character(df[,C_col]) | is.factor(df[,C_col])) {
		gp = gp + scale_color_manual(values=palette)
	} else {
		gp = gp + scale_color_gradientn(colours=palette)
	}
}

opt$x_title = gsub('\\\\n', "\n", opt$x_title)
opt$y_title = gsub('\\\\n', "\n", opt$y_title)
gp = gp + labs(
	x=opt$x_title, 
	y=opt$y_title, 
	title=plot_title
)

# Facet

if (!is.null(opt$facet_by)) {
	facet_col = colnames(df)[opt$facet_by]
	facet_formula = as.formula(sprintf("~%s", facet_col))
	gp = gp + facet_wrap(facet_formula, scale=opt$facet_scale, nrow=opt$facet_nrow)
}


# Legend
gp = gp + guides(colour = guide_legend(override.aes = list(size=5, alpha=1))) 


# Add lines
if (!is.null(opt$abline)) {
	ablines = strsplit(opt$abline, ";")[[1]]
	for (abline in ablines) {
		in_sl = as.numeric(strsplit(abline, ",")[[1]])
		intercept = in_sl[1]
		slope = in_sl[2]
		gp = gp + geom_abline(intercept=intercept, slope=slope, color="grey")
	}
}

# Add the diagonal line
if (opt$diagonal) {
	gp = gp + geom_abline(intercept=0, slope=1, color="grey")
}

# Add the off-diagonal line
if (opt$off_diagonal) {
	gp = gp + geom_abline(intercept=1, slope=-1, color="grey")
}

# Add the regression line
if (opt$linear_regression) {

	if (!is.null(opt$facet_by)) {	
		corDf = data.frame( 
			t( sapply(split(df, df[[facet_col]]), function(d) 
				cor.test(d[[x_col]], d[[y_col]])[c("estimate", "p.value")]
			)
			)
		)
		corDf[[facet_col]] = rownames(corDf)
		corDf[["label"]] = apply(corDf, 1, function(r) {sprintf("c=%.2f\npv=%s", r[1], format(r[2], digits=2))})
	}

	if (opt$verbose) {print(head(df))}
	if (opt$x_log) {
		x_col = sprintf("log10(%s)", x_col)
	}
	if (opt$y_log) {
		y_col = sprintf("log10(%s)", y_col)
	}
	formula = as.formula(sprintf("%s~%s", y_col, x_col))
	if (!is.null(opt$facet_by)) {
		formula = as.formula(sprintf("%s~%s+%s", y_col, x_col, facet_col))
	}
	coeff = lm(formula, df)$coefficients
#	gp = gp + geom_abline(intercept=coeff[1], slope=coeff[2])
	gp = gp + geom_smooth(method = "lm", se=FALSE, color="black", mapping=aes_string(x_col, y_col))

	gp = gp + geom_text(data=corDf, aes(x=min(df[[x_col]]), y=max(df[[y_col]]), label=label), vjust=1, hjust=0)
}






# Change to log scale

if (opt$x_log) {gp = gp + scale_x_log10()}
if (opt$y_log) {gp = gp + scale_y_log10()}

# x-axis limits
x_limits = NULL
if (!is.null(opt$x_limits)) {
    opt$x_limits = gsub("\\", "", opt$x_limits, fixed=TRUE)
    x_limits = as.numeric(strsplit(opt$x_limits, ",")[[1]])
	gp = gp + scale_x_continuous(limits=x_limits)	
}

# y-ayis limits
y_limits = NULL
if (!is.null(opt$y_limits)) {
    opt$y_limits = gsub("\\", "", opt$y_limits, fixed=TRUE)
    y_limits = as.numeric(strsplit(opt$y_limits, ",")[[1]])
	gp = gp + scale_y_continuous(limits=y_limits)	
}




ggsave(opt$output, h=opt$height, w=opt$width)
	




q(save='no')

