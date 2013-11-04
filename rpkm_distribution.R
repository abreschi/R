#!/usr/bin/env Rscript

##------------
## LIBRARIES
##------------ 

cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(plyr))
cat("DONE\n\n")

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=%default]"),
make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment. If empty, column names are used as labels."),
make_option(c("-r", "--representation"), help="choose the representation <boxplot>, <density>, <histogram>, <boxviol> [default=%default]", default="boxplot"),
make_option(c("-o", "--output"), help="additional tags for otuput", default='out'),
make_option(c("-f", "--fill_by"), help="choose the color you want to fill by. Leave empty for no fill", type='character'),
make_option(c("-a", "--alpha_by"), help="choose the fator for the transparency in boxplot. Leave empty for no transparency", type="character"),
make_option(c("-v", "--value"), help="give a name for the value you are plotting the distribution of [default=%default]", default="rpkm"), 
make_option(c("-w", "--wrap"), action="store_true", help="say if you want the density plot in facet_wrap [default=%default]", default=FALSE), 
make_option(c("-T", "--title"), help="give a title to the plot [default=%default]", default=""), 
make_option(c("-P", "--palette"), help="palette name [default=%default]", default="cbbPalette"),
make_option(c("-G", "--no_guide"), help="use this to remove the legend from the plot [default=%default]", default=FALSE, action="store_true"),
make_option(c("-t", "--tags"), help="comma-separated field names you want to display in the labels. Leave default for using column names [default=%default]", default="labExpId")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)

##--------------------##
## CLUSTERING SAMPLES ##
##--------------------##

# read the matrix from the command line
m = read.table(opt$input_matrix, h=T)
ylab = opt$value

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# substitute the matrix with its log if required by the user
if (opt$log) {
m = log10(replace(m, is.na(m), 0) + opt$pseudocount);
ylab = sprintf('log10(%s+%s)', opt$value, opt$pseudocount)}

# prepare data.frame for ggplot
tags = strsplit(opt$tags, ",")[[1]]
df = melt(m, variable.name = "labExpId", value.name=opt$value)

# read the metadata from the metadata file if present
if (!is.null(opt$metadata)) {

mdata = read.table(opt$metadata, h=T, sep='\t')
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))

# prepare data.frame for ggplot
df = merge(unique(mdata[c("labExpId", tags, opt$fill_by, opt$alpha_by)]), df, by="labExpId")

# change color palette in case lenght is too much
if (length(unique(df[,opt$fill_by])) > length(cbbPalette)) {cbbPalette <- rainbow(length(unique(df[,opt$fill_by])))}
}


df$labels = apply(df[tags], 1, paste, collapse="_")
 

###############
# OUTPUT 
###############

output = sprintf("%s.log_%s.psd_%s.%s", opt$representation, 
ifelse(opt$log, "T", "F"), ifelse(opt$log,opt$pseudocount,"NA"),opt$output)

# plotting...

theme_set(theme_bw())

if (opt$representation == "boxplot") {
pdf(sprintf("%s.pdf", output), h=log2(ncol(m)), w=10)
gp = ggplot(df, aes_string(x="labels", y=opt$value))
gp = gp + geom_boxplot(outlier.size=0.5, size=0.5, aes_string(fill=opt$fill_by, alpha=opt$alpha_by))
gp = gp + labs(y=ylab, x='', title=opt$title)
gp = gp + coord_flip()
if(opt$palette == "cbbPalette"){gp = gp + scale_fill_manual(values = cbbPalette)} else {
gp = gp + scale_fill_brewer(palette=opt$palette)}
gp = gp + scale_alpha_manual(values=c(.5,1))
gp = gp + scale_y_continuous(breaks = floor(min(df[is.finite(df[,opt$value]),opt$value])):ceiling(max(df[is.finite(df[,opt$value]),opt$value])))
gp = gp + theme(axis.text=element_text(size=40/log10(length(df$labels))))
}

if (opt$representation == "boxviol") {
pdf(sprintf("%s.pdf", output), h=log2(ncol(m)), w=10)
gp = ggplot(df, aes_string(x="labels", y=opt$value))
gp = gp + geom_violin(size=.4,aes_string(fill=opt$fill_by))
gp = gp + geom_boxplot(alpha=0.2, size=.2)
gp = gp + labs(y=ylab, x='', title=opt$title)
gp = gp + coord_flip()
gp = gp + scale_fill_manual(values = cbbPalette)
gp = gp + stat_summary(fun.y="mean", colour="red", geom="point", size = 3)
gp = gp + scale_y_continuous(breaks = floor(min(df[is.finite(df[,opt$value]),opt$value])):ceiling(max(df[is.finite(df[,opt$value]),opt$value])))
gp = gp + theme(axis.text=element_text(size=40/log10(length(df$labels))))
}

if (opt$representation == "density") {
h=5; w=5
gp = ggplot(df, aes_string(x=opt$value))
if (opt$wrap) {
h=7; w=7
gp = gp + geom_density(aes_string(color=if(is.na(opt$fill_by)){NULL}else{opt$fill_by}))
gp = gp + facet_wrap(~labels)} else {
gp = gp + geom_density( aes_string(group="labels", color=if(is.na(opt$fill_by)){NULL}else{opt$fill_by}) )}
gp = gp + labs(x=ylab, title=opt$title)
gp = gp + scale_color_manual(values = cbbPalette)
gp = gp + scale_x_continuous(breaks = floor(min(df[is.finite(df[,opt$value]),opt$value])):ceiling(max(df[is.finite(df[,opt$value]),opt$value])))
gp = gp + theme(axis.text=element_text(size=40/log10(length(df$labels))))
pdf(sprintf("%s.pdf", output), h=h, w=w)
}

if (opt$representation == "histogram") {
pdf(sprintf("%s.pdf", output))
df$labels = apply(df[tags], 1, paste, collapse="\n") # only in this case collapse with new line
gp = ggplot(df, aes_string(x=opt$value))
gp = gp + geom_histogram(aes_string(fill=if(is.na(opt$fill_by)){NULL}else{opt$fill_by}, y='..count..'), right=TRUE, 
origin=min(df[is.finite(df[,opt$value]),opt$value]))
gp = gp + facet_wrap(~labels)
gp = gp + labs(x=ylab, y='', title=opt$title)
gp = gp + scale_fill_manual(values = cbbPalette)
gp = gp + scale_x_continuous(breaks = floor(min(df[is.finite(df[,opt$value]),opt$value])):ceiling(max(df[is.finite(df[,opt$value]),opt$value])))
gp = gp + theme(axis.text=element_text(size=40/log10(length(df$labels))))
}

if (opt$no_guide) {gp = gp + theme(legend.position="none")}

gp

dev.off()

q(save='no')
