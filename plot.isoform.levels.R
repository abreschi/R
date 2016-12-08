#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))

options(stringsAsFactors=FALSE)


##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix with the transcript rpkms and gene id, columns names are trid, gnid, labExpId1, labExpId2, ..."),
make_option(c("-l", "--log"), help="apply the log"),
make_option(c("-G", "--gene"), help="the gene of interest"),
make_option(c("-o", "--output"), default="isoforms.pdf", help="choose the name for the output file (with extension) [default=%default]"),
make_option(c("-m", "--metadata"), help="metadata file"),
make_option(c("--merge_mdata_on"), default="labExpId", help="metadata column with headers of the input matrix [default=%default]"),
make_option(c("-f", "--field"), default="cell", help="column name of label [default=%default]"),
make_option(c("-c", "--colors"), default="~/R/palettes/rainbow.6.txt", help="color file, the format is RGB [default=%default]"),
make_option(c("-W", "--width"), default=7, help="width of the plot, in inches [default=%default]"),
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

palette = as.character(read.table(opt$colors, h=F, comment.char="%")$V1)

m = read.table(opt$input_matrix, h=T)
colnames(m)[1:2] <- c("trid", "gnid")
m = subset(m, gnid == opt$gene)

if (nrow(m) == 1) {
	cat("The gene has only one transcript! \n")
	cat("EXIT \n")
	q(save="no")
}

m = melt(m)

m = ddply(m, .(gnid, variable), transform, sum=sum(value, na.rm=T))
m$ratio = with(m, value/sum*100)
m$sum = signif(m[,"sum"], 3)

print(head(m))

# Read metadata
if (!is.null(opt$metadata)) {
	field = opt$field
#	opt$merge_mdata_on = "labExpId"
	mdata = read.table(opt$metadata, h=T, sep="\t")
	mdata[,opt$merge_mdata_on] <- gsub(",",".", mdata[,opt$merge_mdata_on])
	mdata = unique(mdata[,unique(c(field, opt$merge_mdata_on))])
	m = merge(m, mdata, by.x="variable", by.y=opt$merge_mdata_on)
}

print(head(m))

if (opt$representation == "stack") {
	gp = ggplot(m, aes_string(x=field, y="ratio")) 
#	gp = gp + geom_bar(position=position_stack(width=1), stat="identity", aes(fill=trid), color="black")
	gp = gp + geom_bar(position="stack", stat="identity", aes(fill=trid), color="black")
	gp = gp + geom_text(aes_string(x=field, y=102, label="sum"), angle=45, hjust=0)
	gp = gp + scale_fill_manual(values=palette, opt$gene)
	gp = gp + theme(axis.text.x = element_text(angle=45, hjust=1))
	gp = gp + labs(x="")
	gp = gp + scale_y_continuous(limits=c(0,105))
}


ggsave(opt$output, w=opt$width, h=6)

q(save="no")
