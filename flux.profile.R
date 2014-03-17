#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input"),
	help="Flux profile, .json"),

make_option(c("-o", "--output"), default="flux.profile.out.pdf",
	help="Output file name [default=%default]"),

make_option(c("-t", "--title"), default="",
	help="Plot title [default=%default]"),

make_option(c("--palette"), type="character", default="/users/rg/abreschi/R/palettes/rainbow.2.txt",
	help="file with the palette [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="if you want more output [default=%default]")

)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

#################
# LIBRARIES 
#################

suppressPackageStartupMessages(library("RJSONIO"))
suppressPackageStartupMessages(library("ggplot2"))



# BEGIN

palette = read.table(opt$palette, h=F, comment.char="%")$V1


m = fromJSON(opt$input)


dfSense = (do.call(
	rbind, sapply(
		1:5, function(i) 
		data.frame(
			x=seq_along(m$masters[[i]]$sense), 
			y=m$masters[[i]]$sense,
			strand="sense",
			bin=i
			), 
		simplify=F
		)
	)
)

dfAntisense = (do.call(
	rbind, sapply(
		1:5, function(i) 
		data.frame(
			x=seq_along(m$masters[[i]]$asense), 
			y=m$masters[[i]]$asense,
			strand="antisense",
			bin=i
			), 
		simplify=F
		)
	)
)


df = rbind(dfSense, dfAntisense)


# PLOT
theme_set(theme_bw(base_size=16))

gp = ggplot(df, aes(x=x, y=y))
gp = gp + geom_line(aes(group=bin, color=strand))
gp = gp + facet_wrap(~bin, scale="free_x", nrow=1)
gp = gp + scale_y_log10()
gp = gp + labs(y="log10(signal)", x="position", title=opt$title)
gp = gp + scale_color_manual(values=palette)
gp = gp + theme(axis.text.x=element_text(angle=45, hjust=1))


ggsave(opt$output, w=10, h=3)

q(save="no")

