#!/usr/bin/env Rscript


# Seed for jittering
set.seed(123)


#~~~~~~~~~~~~
# Libraries
#~~~~~~~~~~~~

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("grid"))


# OPTION PARSING

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="input file. [default=%default]"),

make_option(c("--header"), default=FALSE,
	help="the file has header [deafult=%default]"),

make_option(c("-o", "--output"), default="ggpie.out.pdf",
	help="output file name. Must have a proper image extension (e.g. .pdf, .png) [deafult=%default]"),

make_option(c("-p", "--percentages"), type='integer',
	help="column with the proportions as percentages"),

make_option(c("-n", "--counts"), type='integer',
	help="column with the counts. If you provide counts instead of percentages"),

make_option(c("-f", "--fill_by"), default=2, type="integer",
	help="column with the levels for filling [default=%default]"),

make_option(c("-P", "--palette"), default="/users/rg/abreschi/R/palettes/Set2.8.txt",
	help="palette file [default=%default]"),

make_option(c("-T", "--ggtitle"),
	help="Main title. Leave empty for no title"),

make_option(c("-t", "--fill_title"), 
	help="title for the fill legend")

)


parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

#~~~~~~~~~~~~
# BEGIN
#~~~~~~~~~~~~

# Check that counts and percentages are not provided simultaneously
if (!is.null(opt$percentages) & !is.null(opt$counts)) {
	cat("ERROR: cannot provide both percentages and counts\nAborted\n")
	q(save='no')
}

if (is.null(opt$percentages) & is.null(opt$counts)) {
	cat("MISSING INPUT: provide either percentages or counts\nAborted\n")
	q(save='no')
}
	


# read input
if (opt$input == "stdin") {
	m = read.table(file("stdin"), h=opt$header, sep="\t")
} else {
	m = read.table(opt$input, h=opt$header, sep="\t")
}

# read palette file
palette = as.character(read.table(opt$palette, h=FALSE, comment.char="%")$V1)

# Compute proportions if counts are given
if (!is.null(opt$counts)) {
	counts = colnames(m)[opt$counts]
	m$perc = m[,counts]/sum(m[,counts])*100
	opt$percentages = which(colnames(m) == "perc") 
}


#~~~~~~~~~~~
# Plot
#~~~~~~~~~~~

theme_set(theme_bw(base_size=20))

theme_update(
	panel.grid = element_blank(),
	panel.border = element_blank(),
	axis.ticks = element_blank(),
	axis.ticks.margin = unit(0,'in'),
	axis.title = element_blank(),
	plot.margin = unit(c(0.01, 0.01, 0.01, 0.01),'in'),
	legend.key = element_blank()
)



l = cumsum(m[,opt$percentages]) - c(m[1,opt$percentages], diff(cumsum(m[,opt$percentages])))/2
x_labels = paste(round(m[,opt$percentages], 1), "%", sep="")


gp = ggplot(m, aes_string(x=1, y=colnames(m)[opt$percentages]))
gp = gp + geom_bar(stat="identity", position="stack", aes_string(fill=colnames(m)[opt$fill_by]))
gp = gp + coord_polar(theta="y")
gp = gp + scale_fill_manual(values=palette, name=opt$fill_title)
gp = gp + scale_x_continuous(labels=NULL)
#gp = gp + scale_y_continuous(labels=x_labels, breaks=l)
gp = gp + scale_y_continuous(labels=NULL, breaks=NULL)
gp = gp + geom_text(aes(x=1.5, y=l, label=x_labels), position=position_jitter(w=0.0, h=2.0))
gp = gp + labs(title=opt$ggtitle)


ggsave(opt$output, w=6, h=3.5)

q(save='no')
