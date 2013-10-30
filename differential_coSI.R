
##------------
## LIBRARIES
##------------ 

cat("Loading libraries...")

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
#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
#make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-o", "--output"), help="additional tags for otuput", default="out"),
make_option(c("-c", "--diff_by"), help="choose the factor you want to differ by. Only one factor", type='character')
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)



##--------------------##
## CLUSTERING SAMPLES ##
##--------------------##
output = sprintf("coSI_DE.%s", opt$output)

# 1. read the matrix from the command line
m = read.table(opt$input_matrix, h=T)

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# 2. remove missing values
m <- na.omit(m)

# 3. read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')
mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))
if (!is.null(opt$diff_by)) {opt$diff_by = strsplit(opt$diff_by, ",")[[1]]}

# 4. prepare data.frame for ggplot
df = melt(as.matrix(m), varnames=c("variable", "labExpId"), value.name="value")
df = merge(unique(mdata[unique(c("labExpId", opt$diff_by))]), df, by="labExpId")

# 5. do the average according to factor
new_df = aggregate( as.formula(sprintf("value~%s+variable", opt$diff_by)), df, FUN=mean, na.action=na.omit)
new_df = dcast(new_df, as.formula(sprintf("variable~%s", opt$diff_by)))

# compute M and A
new_df$A = 0.5*(new_df[,2]+new_df[,3])
new_df$M = (new_df[,3]) - (new_df[,2])


# Randomization
set.seed(123)

dfr = df
dfr$value = sample(df$value)
new_dfr = aggregate( as.formula(sprintf("value~%s+variable", opt$diff_by)), dfr, FUN=mean, na.action=na.omit)
new_dfr = dcast(new_dfr, as.formula(sprintf("variable~%s", opt$diff_by)))

new_dfr$A = 0.5*(new_dfr[,2]+new_dfr[,3])
new_dfr$M = (new_dfr[,3])-(new_dfr[,2])


# Add statistical significance after randomizing
new_df$Z = (new_df$M - mean(new_dfr$M))/sd(new_dfr$M)
new_df$pv = 1 - pnorm(abs(new_df$Z))

 

###############
# OUTPUT 
###############

# WRITE TABLE

write.table(new_df, sprintf("%s.tsv", output), row.names=F, quote=F, sep="\t")

# plotting...

pdf(sprintf("%s.pdf",output), h=5, w=6)

theme_set(theme_bw())

gp = ggplot(new_df, aes(x=A, y=M)) 
gp = gp + geom_point(aes(color=cut(abs(Z), breaks =c(0, 1, 1.5, 2, 2.5 ,3 , Inf))))
gp = gp + labs(y="cosi1-cosi2", x='1/2(cosi1+cosi2)')
gp = gp + scale_color_manual(name = "Z-score", values = cbbPalette)
gp = gp + scale_y_continuous(expand=c(0.01,0))
gp = gp + theme(axis.text = element_text(size=15))
gp

gp = ggplot(new_df, aes(x=A, y=M)) 
gp = gp + geom_point(aes(color=cut(abs(pv), breaks =c(-Inf, 0, 0.001, 0.01, 0.05 ,  1))))
gp = gp + labs(y="cosi1-cosi2", x='1/2(cosi1+cosi2)')
gp = gp + scale_color_manual(name = "p-value", values = cbbPalette)
gp = gp + scale_y_continuous(expand=c(0.01,0))
gp = gp + theme(axis.text = element_text(size=15))
gp

dev.off()

q(save='no')
