#!/usr/bin/env Rscript

##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(plyr))


#options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-o", "--output"), help="additional flags for otuput", default="out"),
make_option(c("-W", "--width"), default=5, help="width in inches [default=%default]"),
make_option(c("-H", "--height"), default=5, help="height in inches [default=%default]"),
#make_option(c("-c", "--color_by"), help="choose the color you want to color by [default=NA]", type='character', default=NA),
make_option(c("-f", "--field"), help="dashboard field by which the individuals are grouped")
#make_option(c("-t", "--tags"), help="comma-separated field names you want to display in the labels", default="cell,sex,age")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)

na2null = function(x) if(is.na(x)) {return(NULL)}else{return(x)}


##--------------------##
## CLUSTERING SAMPLES ##
##--------------------##
output = sprintf("%s.%s", basename(opt$input_matrix), opt$output)

# read the matrix from the command line
m = read.table(opt$input_matrix, h=F, col.names=c("element","labExpId","n_det_el","prop","tag"))

# read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))

# prepare data.frame for ggplot
df = merge(subset(m, tag=="individual"), unique(mdata[c("labExpId", opt$field)]), by='labExpId')


# duplicate the data.frame to plot all if an opt$field is specified
if (!is.null(opt$field)) {
df_copy = df;
df_copy[opt$field]= "all";
x_axis = opt$field
# attach the two data frames
new_df = rbind(df, df_copy)
}else{
x_axis='x_lab';
df$x_lab = 'all'
new_df = df}


# order the elements with the hierarchy: junction, exon, transcript, gene
new_df$element <- factor(new_df$element, levels=c('junction','exon','transcript','gene'))

n_colors = length(unique(subset(m, tag=='cumulative')$labExpId))
gp = ggplot(new_df, aes_string(x=x_axis, y="prop"))
gp = gp + geom_boxplot() 
gp = gp + facet_grid(element~.) 
gp = gp + geom_point(data=subset(m, tag=='cumulative'), aes_string(x="labExpId", y="prop"), size = 4, alpha = 0.7, col=cbbPalette[2:2+n_colors])
gp = gp + labs(y='Proportion of detected elements (%)', x="")
gp = gp + theme(axis.text = element_text(size=13))
gp = gp + theme(axis.text.x = element_text(angle=45, hjust=1))
#gp

w=opt$width
h=opt$height

ggsave(filename=sprintf("%s.pdf", output), h=h, w=w)
ggsave(filename=sprintf("%s.png", output), h=h, w=w)
ggsave(filename=sprintf("%s.eps", output), h=h, w=w)


q(save='no')
