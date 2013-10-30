

##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library(reshape))
#suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))

options(stringsAsFactors=F)
pseudocount = 1e-04
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#8470FF", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
custom_palette <- rgb(matrix(c(0,0,0,0,73,73,0,146,146,255,109,182,255,182,119,73,0,146,0,109,219,182,109,255,109,182,255,182,219,255,146,0,0,146,73,0,219,209,0,36,255,36,255,255,109), ncol=3, byrow=T), max=255)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

my_palette <- cbbPalette

##################
# OPTION PARSING
##################

option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze"),
make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-m", "--metadata"), help="A list of tsv files with metadata on matrix experiment.\n\t\tThey must be in the format 'file1.tsv,file2.tsv' and contain a key column named 'labExpId'. Can be omitted"),
make_option(c("-o", "--output"), help="additional info you want to put in the output file name", default="out"),
make_option(c("-c", "--color_by"), help="choose the fields in the metadata you want to color by", type='character', default=NA),
make_option(c("-s", "--shape_by"), help="choose the fields in the metadata you want to shape by", type='character', default=NA),
make_option(c("-r", "--row_as_variables"), action="store_true", help="select this if you want rows as variables [default=%default]", default=FALSE),
make_option(c("-C", "--princomp"), help="choose the principal components you want to plot. With 3 PC it gives a 3d plot [default='PC1,PC2']", default="PC1,PC2")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

print(opt)


###############
# BEGIN
##############

# read input tables
m = read.table(opt$input_matrix, h=T)


# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# apply the log if required
if (opt$log) {m = log2(replace(m, m==0 | is.na(m), pseudocount))}

# apply pca
if (opt$row_as_variable) {
m_pca = prcomp(na.omit(m), center=FALSE, scale.=FALSE)}else{
m_pca = prcomp(t(na.omit(m)), center=FALSE, scale.=FALSE)}

print(dim(na.omit(m)))

#
# HANDLING METADATA
#
# add metadata to pca results, they should be input in the command line in the future
if (is.na(opt$color_by)) {color_by=NULL
}else{color_by = color_by = strsplit(opt$color_by, ",")[[1]]}
if (is.na(opt$shape_by)) {shape_by=NULL
}else{shape_by = strsplit(opt$shape_by, ",")[[1]]}

# read metadata, one or more table to be merged on labExpId
if (!is.null(opt$metadata)){
	metadata = strsplit(opt$metadata, ",")[[1]]
	for (i in seq_along(metadata)) {
		mdata = read.table(metadata[i], h=T, sep="\t", row.names=NULL);
		if ('labExpId' %in% colnames(mdata)) {
		mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))}
		if ( i==1 ) {
			new_mdata = mdata
		}else{
		new_mdata = merge(mdata, new_mdata, by=c('cell'))
		}
	}

cat('append metadata...')

df = merge(as.data.frame(m_pca$x),
unique(new_mdata[c(color_by, shape_by, 'labExpId')]), by.x='row.names', by.y='labExpId', all.x=T)
}else{
df = as.data.frame(m_pca$x)}



#dfl = as.data.frame(m_pca$rotation)


#########
# OUTPUT
#########

output_name = sprintf("log_%s.pseudo_%s.colby_%s.shpby_%s.prcomp_%s.%s.pca", 
opt$log, opt$pseudocount, gsub(",",".",opt$color_by), gsub(",", ".", opt$shape_by), gsub(",",".",opt$princomp), opt$output)

write.table(m_pca$x, sprintf("%s.tsv", output_name), quote=F, sep="\t")

prinComp = strsplit(opt$princomp, ",")[[1]]
variance_x = round(summary(m_pca)$importance[2,as.numeric(gsub("PC","",prinComp[1]))],3)*100
variance_y = round(summary(m_pca)$importance[2,as.numeric(gsub("PC","",prinComp[2]))],3)*100

pdf(sprintf("%s.pdf", output_name), w=9)
# plot parameters
pts = 5
shapes = c(16, 15, 0, 14:1, 18, 17)


theme_set(theme_bw())

if (length(prinComp) == 2){
# plotting...
gp = ggplot(df, aes_string(x=prinComp[1],y=prinComp[2]));
if (!is.na(opt$color_by)) {gp_color_by=interaction(df[color_by])} else {gp_color_by=NULL}
if (!is.na(opt$shape_by)) {gp_shape_by=interaction(df[shape_by]);
gp_shape_by <- factor(gp_shape_by, levels=sort(levels(gp_shape_by)))} else {gp_shape_by=NULL}
gp = gp + geom_point(aes(col=gp_color_by, shape=gp_shape_by), size=pts);
gp = gp + labs(title="");
gp = gp + labs(x=sprintf('%s (%s%%)', prinComp[1], variance_x));
gp = gp + labs(y=sprintf('%s (%s%%)', prinComp[2],variance_y));
gp = gp + scale_color_manual(name=opt$color_by, values=my_palette)
gp = gp + scale_shape_manual(name='Sample', values=shapes);
gp} 




# ---------------------------------------
# 3d scatterplot

if (length(prinComp) == 3) {

par(xpd=NA, omi=c(.5,.5,.5,.5))

variance_z = round(summary(m_pca)$importance[2,as.numeric(gsub("PC","",prinComp[3]))],3)*100
suppressPackageStartupMessages(library(scatterplot3d))
if (!is.na(opt$color_by)) {gp_color=my_palette[interaction(df[color_by])]} else {gp_color="black"}
if (!is.na(opt$shape_by)) {gp_shape_by=interaction(df[shape_by]);
gp_shape_by <- factor(gp_shape_by, levels=sort(intersect(levels(gp_shape_by), gp_shape_by))); gp_shape=shapes[gp_shape_by]} else {gp_shape_by=NULL}
scatterplot3d(df[prinComp], 
color = gp_color,
pch = gp_shape,
xlab = sprintf('%s (%s%%)', prinComp[1], variance_x),
ylab = sprintf('%s (%s%%)', prinComp[2], variance_y),
zlab = sprintf('%s (%s%%)', prinComp[3], variance_z)
)

if (!is.na(opt$color_by)) {
legend(log(max(df[prinComp[1]])),7.2,levels(interaction(df[color_by])), 
fill=my_palette[seq_along(levels(interaction(df[color_by])))])}

if (!is.na(opt$shape_by)) {
legend(-log(abs(min(df[prinComp[1]])))+1.5,7.2,levels(gp_shape_by), 
pch=shapes[seq_along(levels(gp_shape_by))])}
}


dev.off()
q(save='no')

