#!/usr/bin/env Rscript

##------------
## LIBRARIES
##------------ 
cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
#suppressPackageStartupMessages(library(plyr))
cat("DONE\n\n")

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################

option_list <- list(

make_option(c("-i", "--input_files"), type="character", 
	help="list of files with the bigWig profiles, comma-separated. If only one file is given, can be stdin"),
make_option(c("--header"), action="store_true", default=FALSE, 
	help="The files have header [default=%default]"),
make_option(c("-L", "--labels"), help="list of labels with the labels of each file, commma-separated.\n
They must be in the same order as the file list", type="character"),
make_option(c("-v", "--vertical_lines"), help="specify where you want the vertical lines [default=%default]", type="character", default="0"),
#make_option(c("-f", "--facet"), type="integer", help="column index to facet"),
make_option(c("-t", "--title"), help="Main title for the plot [default=%default]", default=""),
make_option(c("-o", "--output"), help="output file name with extension [default=%default]", default="profile.pdf"),
make_option(c("-y", "--y_title"), help="title for the y-axis [default=%default]", default="norm_read_density")
#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
#print(opt)

# ========
# BEGIN
# ========
if (!is.null(opt$labels)) {labels = strsplit(opt$labels, ",")[[1]]}

vertical_lines = as.numeric(strsplit(opt$vertical_lines, ",")[[1]])


# Read data

if (opt$input_files != "stdin") {
	input_files = strsplit(opt$input_files, ",")[[1]]
	i = 1
	for (f in input_files) {
		if (i == 1 ) {
			m = read.table(f, h=opt$header)
			colnames(m)[1] <- "position"
			colnames(m)[2] <- labels[i]
		}
		if (i != 1) {
			tmp = read.table(f, h=F)
			colnames(tmp)[1] <- "position"
			colnames(tmp)[2] <- labels[i]
			m = merge(m, tmp, by = "position")
		}
		i = i+1
	}
} else {
	m = read.table(file('stdin'), h=opt$header)
	colnames(m)[1] <- "position"
	if (is.null(opt$labels)) {
		colnames(m)[2] <- 'label'
	} else {
		colnames(m)[2] <- labels[1]
	}
}


# Melt data

df = melt(m, id.vars="position")


# GGPLOT

theme_set(theme_bw(base_size=20))

gp = ggplot(df, aes(x=position, y=value)) 
gp = gp + geom_line(aes(group=variable, color=variable), alpha=0.5, size=1)
gp = gp + geom_vline(xintercept=vertical_lines, color="grey", linetype="longdash")
gp = gp + labs(title=opt$title, y=opt$y_title)
gp = gp + scale_color_manual(values=cbbPalette)
if (is.null(opt$label)) {gp = gp + theme(legend.position='none')}

ggsave(opt$output, h=5, w=7)

# EXIT
quit(save='no')
