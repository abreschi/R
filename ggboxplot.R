#!/usr/bin/env Rscript

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"), type="character", default='stdin',
	help="tab-separated file. Can be stdin [default=%default]"),

make_option(c("-o", "--output"), default="ggboxplot.pdf", 
	help="output file name with extension [default=%default]"),

make_option(c("--header"), action="store_true", default=FALSE, 
	help="The file has header [default=%default]"),

make_option(c("-c", "--color_by"), type='numeric',
	help="Column index with the color factor. Leave empty for no color"),

make_option(c("-f", "--fill_by"), type='numeric', 
	help="Column index of the fill factor, Leave empty for no fill"),

make_option(c("--xy"), type='character', default="1,2",
	help="Column indeces of the x and y axes, respectively, comma-separated [default=%default]"),

make_option(c("-r", "--representation"), default="boxplot",
	help="Choose representation [default=%default]"),

make_option(c("--facet_by"), type="integer", 
	help="column index to facet"),

make_option(c("--log"), action="store_true", default=FALSE,
	help="apply the log to the y-axis [default=%default]"),

make_option(c("-p", "--pseudocount"), type="numeric", default=1e-03,
	help="Pseudocount for the log on the y-axis [default=%default]"),

make_option(c("-t", "--title"), type="character", 
	help="Main title for the plot. Leave emtpy for no title"),

make_option(c("--y_title"),  
	help="title for the y-axis. Leave empty for no title"),

make_option(c("--x_title"), 
	help="title for the x-axis. Leave empty for no title"),

make_option(c("-P", "--palette"), 
	help='File with colors for the lines. Leave empty to use even color spacing'),

make_option(c("-v", "--verbose"), action='store_true', default=FALSE,
	help="Verbose output [default=%default]")

#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list = option_list,
	description = "From a column file, plot a column vs another as lines"
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
#suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


# ========
# BEGIN
# ========


# Read data
if (opt$input == "stdin") {input=file('stdin')} else {input=opt$input}
m = read.table(input, h=opt$header)
if(opt$verbose) {print(head(m))}

# Read the axes
xy = strsplit(opt$xy, ",")[[1]]
x_col = as.numeric(xy[1])
y_col = as.numeric(xy[2])
x = colnames(m)[x_col]
y = colnames(m)[y_col]

# Convert to log if asked
if (opt$log) {
	m[,y] <- log10(m[,y]+opt$pseudocount)
}

# Read palette
if (!is.null(opt$palette)) {
	palette = read.table(opt$palette, h=FALSE, comment.char='%')$V1
}


# # Significance
# levs = levels(as.factor(m[,x]))
# pv = sapply(1:(length(levs)-1), function(i) {
# 	distr1 = m[which(m[,x] == levs[i]), y];
# 	distr2 = m[which(m[,x] == levs[i+1]), y];
# 	return(suppressWarnings(ks.test(distr1, distr2)$p.value))
# 	}
# )
# sign_df = data.frame(x=levs[-length(levs)], text=format(pv, digits=2))

#~~~~~~~~~~~~
# GGPLOT
#~~~~~~~~~~~~

theme_set(theme_bw(base_size=20))
theme_update(
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	legend.key = element_blank(),
	plot.title = element_text(vjust=1)
)

if (!is.null(opt$color_by)) {color_by = colnames(m)[opt$color_by]} else {color_by=NULL}
if (!is.null(opt$fill_by)) {fill_by = colnames(m)[opt$fill_by]} else {fill_by=NULL}
alpha = 0.9


# plot 
gp = ggplot(m, aes_string(x=x, y=y))
if (opt$representation == "boxplot") { 
	gp = gp + geom_boxplot(aes_string(color=color_by, fill=fill_by), alpha=alpha, size=1)
}
if (opt$representation == "violin") {
	gp = gp + geom_violin(aes_string(color=color_by, fill=fill_by), alpha=alpha, size=1)
}
gp = gp + labs(title=opt$title, y=opt$y_title, x=opt$x_title)
gp = gp + theme(axis.text.x=element_text(angle=45, hjust=1))

# Color scale
if (!is.null(opt$palette)) {
	gp = gp + scale_fill_manual(values=palette)
} else {
	gp = gp + scale_fill_hue()
}

# Read facet
if (!is.null(opt$facet_by)) {
	facet_col = colnames(m)[opt$facet_by]
	facet_form = as.formula(sprintf("~%s", facet_col))
	gp = gp + facet_wrap(facet_form)
}
#gp = gp + geom_text(data=sign_df, aes(x=x, y=1, label=text), hjust=-0.5)

ggsave(opt$output, h=5, w=7)

# EXIT
quit(save='no')
