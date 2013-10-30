
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
theme_set(theme_bw(base_size=16))


##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-o", "--output"), help="additional flags for otuput", default="out"),
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
output = sprintf("%s.%s.%s", basename(opt$input_matrix), opt$field, opt$output)

# read the matrix from the command line
m = read.table(opt$input_matrix, h=F, col.names=c("labExpId","n_det_el","element"))

# read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))

# separate total from the rest
mtot = subset(m, element=='total')
m = subset(m , element!='total')

# attach the total as additional column
df = merge(m, mtot, by='labExpId')

# attach the metadata
df = merge(df, unique(mdata[c('labExpId', opt$field)]), by='labExpId')
df$element.x <- factor(df$element.x, levels=c('exonic', 'intronic', 'exonic-intronic', 'genic', 'intergenic'))


gp = ggplot(df) 
gp = gp + geom_boxplot(aes(y=n_det_el.x/n_det_el.y*100, x=element.x)) 
gp = gp + facet_grid(~cell)
gp = gp + labs(y='Proportion of mappings (%)', x="")
gp = gp + theme(axis.text = element_text(size=13, angle=45, h=1))

w=5
h=5

ggsave(filename=sprintf("%s.pdf", output), h=h, w=w)
ggsave(filename=sprintf("%s.png", output), h=h, w=w)
ggsave(filename=sprintf("%s.eps", output), h=h, w=w)


q(save='no')
