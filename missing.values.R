
##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(plyr))


options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-o", "--output"), help="additional tags for otuput [default=%default]", default='out'),
make_option(c("-t", "--tags"), help="comma-separated field names you want to display in the labels", default="cell,sex,age")
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

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# substitute the matrix with its log if required by the user
if (opt$log) {m = log2(replace(m, is.na(m), 0) + opt$pseudocount)}

# read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))

# plot the number of missing values and zeros
stat = data.frame(NAs = apply(m,2,function(x) sum(is.na(x))), 
zeros = apply(m,2,function(x) sum(x==0,na.rm=T)), pos_values = apply(m,2,function(x) sum(!is.na(x)&x!=0)))
labels = strsplit(opt$tags, ',')[[1]]
df = merge(unique(mdata[c('labExpId', labels)]), stat, by.y='row.names', by.x='labExpId')
df = melt(df, ids = c('labExpId',labels), variable = 'values', value.name="Number_of_genes")
df$labels = apply(df[labels], 1, paste, collapse='_' )


output = sprintf("log_%s.psd_%s.%s.missing.values", opt$log, opt$pseudocount, opt$output)
pdf(sprintf('%s.pdf', output), h=log2(ncol(m)), w=7)
gp = ggplot(df, aes(x=labExpId, y=Number_of_genes)) 
gp = gp + geom_bar(aes(fill=values), stat='identity') 
gp = gp + coord_flip()
gp = gp + scale_x_discrete(labels=df$labels)
gp = gp + labs(x='')
gp = gp + theme(axis.text=element_text(size=10/log10(ncol(m))))
gp
dev.off()

q(save='no')
