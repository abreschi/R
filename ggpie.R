#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))


# OPTION PARSING

option_list <- list(

make_option(c("-i", "--input"), default="stdin",
	help="input file. [default=%default]"),

make_option(c("--header"), default=FALSE,
	help="the file has header [deafult=%default]"),

make_option(c("-o", "--output"), default="out.pdf",
	help="output file name [deafult=%default]"),

make_option(c("-p", "--proportions"), default=1, type="integer",
	help="column with the proportions [default=%default]"),

make_option(c("-f", "--fill_by"), default=2, type="integer",
	help="column with the levels for filling [default=%default]"),

make_option(c("-t", "--fill_title"), 
	help="title for the fill legend")

)


parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options



if (opt$input == "stdin") {
	m = read.table(file("stdin"), h=opt$header)
} else {
	m = read.table(opt$input, h=opt$header)
}


theme_set(theme_bw(base_size=16))

l = cumsum(m[,opt$proportions]) - c(m[1,opt$proportions], diff(cumsum(m[,opt$proportions])))/2

gp = ggplot(m, aes_string(x=1, y=colnames(m)[opt$proportions]))
gp = gp + geom_bar(stat="identity", position="stack", aes_string(fill=colnames(m)[opt$fill_by]))
gp = gp + coord_polar(theta="y")
gp = gp + scale_fill_brewer(palette="Set2", opt$fill_title)
gp = gp + labs(x="", y="")
gp = gp + scale_x_continuous(labels=NULL)
gp = gp + scale_y_continuous(labels=round(m[,opt$proportions], 1), breaks=l)
gp = gp + theme(panel.grid=element_blank())

ggsave(opt$output, w=7, h=6)

