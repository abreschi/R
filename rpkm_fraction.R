
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
#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
#make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-o", "--output"), help="additional tags for otuput", default="out"),
make_option(c("-c", "--color_by"), help="choose the color you want to color by. Leave empty for no color", type='character'),
make_option(c("-y", "--linetype_by"), help="choose the factor you want the linetype by. Leave empty for no linetype", type="character"),
make_option(c("-t", "--tags"), help="choose the factor by which grouping the lines [default=%default]", default="labExpId")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)



##--------------------##
## CLUSTERING SAMPLES ##
##--------------------##
output = sprintf("rpkm_fraction.%s", opt$output)

# 1. read the matrix from the command line
m = read.table(opt$input_matrix, h=T)

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}


# 2. sort all the values for each column
sortm <- apply(m, 2, sort, na.last=T, d=T)
#sortm <- sortm[1:min(which(sortm==0))-1,]

# 3. calculate cumulative sum
cumm <- as.data.frame(apply(sortm, 2, function(x) cumsum(x/sum(x,na.rm=T))))


# 4. read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))
if (!is.null(opt$color_by)) {opt$color_by <- strsplit(opt$color_by, ",")[[1]]}

# prepare data.frame for ggplot
df = melt(cumm, variable.name = "labExpId", value.name="rpkm_fraction")
df = merge(unique(mdata[unique(c("labExpId", strsplit(opt$tags, ",")[[1]], opt$color_by, opt$linetype_by))]), df, by="labExpId")
df$labels = apply(df[strsplit(opt$tags, ",")[[1]]], 1, paste, collapse="_")
# add a column with the x index
df  = ddply(df, .(labels), transform, x=seq_along(labels), y=sort(rpkm_fraction, na.last=T, d=F))


###############
# OUTPUT 
###############

# plotting...
base_size=16
#legend_nrow = 18
legend_nrow = 4
theme_set(theme_bw(base_size=base_size))
legend_text_inch = theme_get()$legend.text$size * base_size / 72.72
add_w = legend_text_inch * max(nchar(df$labels)) * ceiling(length(levels(as.factor(df$labels)))/legend_nrow)

pdf(sprintf("%s.pdf",output), h=5, w=6+add_w, title=output)


gp = ggplot(df, aes(x=x, y=y, group=labels))
if (!is.null(opt$color_by)) {gp_color_by=interaction(df[opt$color_by])} else {gp_color_by=NULL}
if (!is.null(opt$linetype_by)) {gp_linetype_by=interaction(df[opt$linetype_by])} else {gp_linetype_by=NULL}
gp = gp + geom_line(aes(color=gp_color_by, linetype=gp_linetype_by))
gp = gp + labs(y="Fraction of gene rpkm", x='Number of genes')
gp = gp + scale_color_manual(values = cbbPalette)
gp = gp + scale_linetype_manual(values=c(2,1))
#gp = gp + scale_color_hue(name=paste(opt$color_by, collapse="."))
#gp = gp + guides(col = guide_legend(nrow = legend_nrow))
gp = gp + scale_x_log10(expand=c(0,0))
gp = gp + scale_y_continuous(expand=c(0.01,0))
gp = gp + annotation_logticks(sides="b")
gp

dev.off()

q(save='no')
