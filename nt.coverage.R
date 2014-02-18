
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
make_option(c("-d", "--outdir"), help="directory for the output", default="./"),
#make_option(c("-c", "--color_by"), help="choose the color you want to color by [default=NA]", type='character', default=NA),
make_option(c("-f", "--field"), help="dashboard field by which the individuals are grouped"),
make_option(c("-H", "--height"), default=5, help="output height in inches [default=%default]"),
make_option(c("-W", "--width"), default=9, help="output width in inches [default=%default]"),
make_option(c("--facet_nrow"), default=1, help="number of rows when faceting [default=%default]")

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
m = read.table(opt$input_matrix, h=F, col.names=c("element","labExpId","n_det_el","prop","tag","det"))

# read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))

# prepare data.frame for ggplot
df = merge(subset(m, tag=="individual"), unique(mdata[c("labExpId", opt$field)]), by='labExpId')

# duplicate the data.frame to plot all
df_copy = df
df_copy[opt$field]= "all"

# attach the two data frames
new_df = rbind(df, df_copy)


# order the elements with the hierarchy: exonic, intronic, intergenic, total
new_df$element <- factor(new_df$element, levels=c('exonic','intronic','intergenic','total'))
#new_df$labels <- factor(new_df$det, labels = c("Genome covered", "Covered by RNA-seq", "RNA-seq distribution"))
new_df$facet2 = interaction(new_df$element, new_df$det)

# create a separate data.frame for cumulative coverage
cumul_df = subset(m, tag=='cumulative')
cumul_df$facet2 = interaction(cumul_df$element, cumul_df$det)
cumul_df[opt$field] = cumul_df$labExpId
#cumul_df$cell = cumul_df$labExpId


#gp = ggplot(new_df, aes(x=element, y=prop))
#gp = gp + geom_boxplot() 
#gp = gp + facet_grid(cell~det)
##gp = gp + facet_grid(sprintf(".~%s", opt$field), scales = 'free_y') 
##gp = gp + facet_grid(new_df$cell~facet, scales = 'free_y') 
##gp = gp + geom_point(data=subset(m, tag=='cumulative'), aes_string(x="element", y="prop"), size = 4, alpha = 0.7)
#gp = gp + labs(y='Proportion of detected nucleotides', x="")
#gp = gp + theme(axis.text = element_text(size=15), axis.text.x=element_text(angle=45))
#gp
#
#gp = ggplot(new_df, aes(x=element, y=prop))
#gp = gp + geom_boxplot(aes(col=cell)) 
#gp = gp + facet_grid(.~labels, scales='free_x')
##gp = gp + facet_grid(sprintf(".~%s", opt$field), scales = 'free_y') 
##gp = gp + facet_grid(new_df$cell~facet, scales = 'free_y') 
##gp = gp + geom_point(data=subset(m, tag=='cumulative'), aes_string(x="element", y="prop"), size = 4, alpha = 0.7)
#gp = gp + labs(y='Proportion of detected nucleotides', x="Genomic domain (D)")
#gp = gp + theme(axis.text = element_text(size=15), axis.text.x=element_text(angle=45, hjust=1))
#gp
#
## almost final without cumulative
#gp = ggplot(new_df, aes(x=facet2, y=prop))
#gp = gp + geom_boxplot(aes(col=labels)) 
#gp = gp + facet_grid(.~cell, scales='free_x')
##gp = gp + facet_grid(sprintf(".~%s", opt$field), scales = 'free_y') 
##gp = gp + facet_grid(new_df$cell~facet, scales = 'free_y') 
##gp = gp + geom_point(data=subset(m, tag=='cumulative'), aes_string(x="element", y="prop"), size = 4, alpha = 0.7)
#gp = gp + labs(y='Proportion of nucleotides (%)', x="Genomic domain (D)")
#gp = gp + theme(axis.text = element_text(size=15), axis.text.x=element_text(angle=45, hjust=1))
#gp = gp + scale_x_discrete(labels=c('total','exonic','intronic','intergenic','exonic','intronic','intergenic'))
#gp

# 

gp = ggplot(new_df, aes(x=facet2, y=prop))
gp = gp + geom_boxplot(aes(fill=as.character(det)) )
#gp = gp + facet_grid(.~cell, scales='free_x')
if (!is.null(opt$field)) {
#gp = gp + facet_grid(sprintf(".~%s", opt$field), scales='free_x')}
gp = gp + facet_wrap(as.formula(sprintf("~%s", opt$field)), nrow=opt$nrow, scales='free_x')}
gp = gp + geom_point(data=cumul_df, aes(x=facet2, y=prop, col=as.character(det)), size = 4, alpha = 0.7)
gp = gp + geom_point(data=cumul_df, aes(x=facet2, y=prop), col='black', size = 4, alpha = 0.7, shape=1)
gp = gp + labs(y='Proportion of nucleotides (%)', x="Genomic domain (D)")
gp = gp + theme(axis.text = element_text(size=15), axis.text.x=element_text(angle=45, hjust=1))
gp = gp + scale_x_discrete(labels=c('total','exonic','intronic','intergenic','exonic','intronic','intergenic'))
gp = gp + scale_fill_hue(name='',labels = c("Genome covered", "Covered by RNA-seq", "RNA-seq distribution"))
gp = gp + scale_color_hue(name='Cumulative',labels = c("Genome covered", "Covered by RNA-seq", "RNA-seq distribution"))


w = opt$width
h = opt$height

ggsave(filename=sprintf("%s/%s.pdf", opt$outdir, output), w=w, h=h)
ggsave(filename=sprintf("%s/%s.png", opt$outdir, output), w=w, h=h)
ggsave(filename=sprintf("%s/%s.eps", opt$outdir, output), w=w, h=h)


q(save='no')
