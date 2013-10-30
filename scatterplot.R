
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
x_psd = 1e-03
y_psd = 1e-03
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), help="the matrix you want to analyze WITH header"),
#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
#make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
make_option(c("-x", "--x_axis"), type='integer', help="the index (1-based) of the column you want on the x axis"),
make_option(c("-y", "--y_axis"), type='integer', help="the index (1-based) of the column you want on the y axis"),
make_option(c("-o", "--output_suffix"), help="additional output strings [default=%default]", default='out'),
make_option(c("-t", "--type"), help="<tile>, <hex>, <scatter> [default=%default]", default="tile"),
make_option(c("-b", "--binwidth"), help="comma-separated values for binwidth x,y [default=%default]", default="1,1"),
make_option(c("--x_log"), action="store_true", help="x values log10 transformed [default=%default]", default=FALSE),
make_option(c("--y_log"), action="store_true", help="y values log10 transformed [default=%default]", default=FALSE),
make_option(c("--x_psd"), help="pseudocount for x values [default=%default]", default=x_psd, type='double'),
make_option(c("--y_psd"), help="pseudocount for y values [default=%default]", default=y_psd, type='double'),
make_option("--x_title", help="write a title for x axis [default=%default]", default="x_title"),
make_option("--y_title", help="write a title for y axis [default=%default]", default="y_title"),
make_option("--legend_title", help="write a title for the legend [default=%default]", default="count"),
make_option("--title", help="write a title for the plot [default=%default]", default="plot_title")
)



parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)



###################
#| BEGIN         |#
###################


m = read.table(opt$input_matrix, h=T)

if (opt$x_log) {m[,opt$x_axis] <- m[,opt$x_axis] + opt$x_psd}
if (opt$y_log) {m[,opt$y_axis] <- m[,opt$y_axis] + opt$y_psd}

df = m

# Pearson correlation coefficient
#corr = round(cor(sapply(df[,opt$x_axis], function(x) ifelse(opt$x_log, log10(x), x)), 
#sapply(df[,opt$y_axis], function(x) ifelse(opt$y_log, log10(x), x)), method='p', use='p'), 2)
pearson = round(cor(sapply(df[,opt$x_axis], function(x) ifelse(opt$x_log, log10(x), x)), 
sapply(df[,opt$y_axis], function(x) ifelse(opt$y_log, log10(x), x)), method='p', use='p'), 2)
spearman = round(cor(sapply(df[,opt$x_axis], function(x) ifelse(opt$x_log, log10(x), x)), 
sapply(df[,opt$y_axis], function(x) ifelse(opt$y_log, log10(x), x)), method='s', use='p'), 2)

# PLOTTING ...

output = sprintf("scatterplot.%s", opt$output_suffix)

theme_set(theme_bw(base_size=16))

bwidth = as.numeric(strsplit(opt$binwidth, ",")[[1]])
x_title = 'Mouse (number of annotated isoforms per gene locus)'
y_title = 'Human (number of annotated isoforms per gene locus)'
plot_title = "Number of annotated isoforms per gene locus"
plot_title = sprintf("%s (p_r=%s; s_r=%s)", opt$title, pearson, spearman)

#countBins = c(0,1,5,10,100,500,1000,2000,3000,Inf)

countBins <- c(0,1,2,5,10,25,50,75,100,500,Inf)

pdf(sprintf("%s.pdf", output))
gp = ggplot(df, aes_string(x=colnames(df[opt$x_axis]), y=colnames(df[opt$y_axis]))) 
if (opt$type == 'tile') {
#gp = gp + stat_bin2d(binwidth=bwidth, aes(fill=cut(..count.., c(0,1,2,5,10,25,50,75,100,500,Inf))),colour="black",size=.2)
gp = gp + stat_bin2d(bins=100)
gp = gp + scale_fill_gradientn(colours=terrain.colors(20), name=opt$legend_title) }
#gp = gp + scale_fill_brewer("count") }
if (opt$type == 'hex') {
gp = gp + geom_hex(aes(fill=cut(..count.., c(0,1,2,5,10,25,50,75,100,500,Inf))), binwidth=bwidth)
gp = gp + scale_fill_manual('counts', values=terrain.colors(length(countBins))) }
if (opt$type == "scatter") {
gp = gp + geom_point(shape=".")
}

#gp = gp + geom_abline(intercept=0, slope=1, color='blue', linetype=4) 
gp = gp + labs(x=opt$x_title, y=opt$y_title, title=plot_title)
#gp = gp + geom_text(aes(x=min(df[opt$x_axis], na.rm=T)+.1, y=max(df[opt$y_axis], na.rm=T)-.1, 
#label=sprintf("p_r=%s\ns_r=%s",pearson, spearman)), color='blue', size=5)

#if (opt$x_log) {gp = gp + scale_x_log10() + annotation_logticks(sides="b")}
#if (opt$y_log) {gp = gp + scale_y_log10() + annotation_logticks(sides="l")}

if (opt$x_log) {gp = gp + scale_x_log10()}
if (opt$y_log) {gp = gp + scale_y_log10()}

gp
dev.off()




q(save='no')

