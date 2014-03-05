#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))

options(stringsAsFactors=FALSE)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix with the transcript rpkms and gene id, columns names are trid, gnid, labExpId1, labExpId2, ..."),
make_option(c("-l", "--log"), help="apply the log"),
make_option(c("-G", "--gene"), help="the gene of interest"),
make_option(c("-o", "--output"), help="choose the name for the output file"),
make_option(c("-c", "--colors"), default="~/R/palettes/rainbow.6.txt", help="color file, the format is RGB [default=%default]"),
make_option(c("-r", "--representation"), default="stack", help="visualization method: stack | [defuault=%default]")

#make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
#make_option(c("-d", "--de"), help='output of edgeR for differentially expressed genes')
)


parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

#print(opt)

#--------------------
# FUNCTIONS
#--------------------

palette = as.character(read.table(opt$colors, h=F)$V1)

m = read.table(opt$input_matrix, h=T)
m = subset(m, gnid == opt$gene)

if (nrow(m) == 1) {
	cat("The gene has only one transcript! \n")
	cat("EXIT \n")
	q(save="no")
}

m = melt(m)

m = ddply(m, .(gnid, variable), transform, sum=sum(value, na.rm=T))
m$ratio = with(m, value/sum*100)



if (opt$representation == "stack") {
	gp = ggplot(m, aes(x=variable, y=ratio)) 
	gp = gp + geom_bar(position=position_stack(width=1), stat="identity", aes(fill=trid), color="black")
	gp = gp + geom_text(aes(x=variable, y=102, label=signif(sum, 3)))
	gp = gp + scale_fill_manual(values=palette, opt$gene)
	gp = gp + theme(axis.text.x = element_text(angle=45, hjust=1))
	gp = gp + labs(x="")
}


ggsave(opt$output, w=7, h=6)

q(save="no")
