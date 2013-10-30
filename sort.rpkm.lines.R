
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
make_option(c("-o", "--output"), help="additional tags for otuput", default="out"),
make_option(c("-c", "--color_by"), help="choose the factor you want to color by", type='character'),
make_option(c("-L", "--linetype_by"), help="choose the factor you want to change the linetype by. Leave empty for none."),
make_option(c("-v", "--value"), help="the value the you are plotting [default=%default]", default="RPKM"),
make_option(c("-t", "--tags"), help="comma-separated field names you want to display in the labels", default="cell,sex,age")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)

na2null = function(x) if(is.na(x)) {return(NULL)}else{return(x)}


##--------------------##
## CLUSTERING SAMPLES ##
##--------------------##
output = sprintf("log_%s.pseudo_%s.colby_%s.%s", opt$log, opt$pseudocount, opt$color_by, opt$output)

# read the matrix from the command line
m = read.table(opt$input_matrix, h=T)
ylab = 'rpkm'

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

#substitute the matrix with its log if required by the user
if (opt$log) {
m = log10(replace(m, is.na(m), 0) + opt$pseudocount);
ylab = sprintf('log10(%s+%s)', opt$value, opt$pseudocount)}

# read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))

# prepare data.frame for ggplot
df = melt(m, variable.name = "labExpId", value.name="rpkm")
df = merge(unique(mdata[unique(c("labExpId", strsplit(opt$tags, ",")[[1]], opt$color_by, opt$linetype_by))]), df, by="labExpId")
df$labels = apply(df[strsplit(opt$tags, ",")[[1]]], 1, paste, collapse="_")
# add a column with the x index
df  = ddply(df, .(labels), transform, x=seq_along(labels), y=sort(rpkm, na.last=T,d=T))
 

###############
# OUTPUT 
###############

# plotting...
pdf(sprintf("%s.pdf",output), h=5, w=5)

theme_set(theme_bw())

gp = ggplot(df, aes(x=x, y=y, group=labels))
gp = gp + geom_line(aes_string(color=opt$color_by, linetype=opt$linetype_by))
gp = gp + labs(y=ylab, x=sprintf('rank(%s)', opt$value))
gp = gp + scale_color_manual(values = cbbPalette)
gp = gp + scale_x_log10(expand=c(0,0))
gp = gp + annotation_logticks(sides="b")
gp = gp + theme(axis.text = element_text(size=15))
gp

dev.off()

q(save='no')
