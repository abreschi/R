#!/usr/bin/env Rscript

# -- Variables --

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-u", "--universe"), help="a list of human gene identifiers (ensEMBL ids), NO header"),
make_option(c("-G", "--genes"), help="a list of human gene identifiers for the foreground (ensEMBL ids), WITH header"),
make_option(c("-c", "--categ"), help="choose the GO category <BP>, <MF>, <CC> [default=%default]", default="BP"),
make_option(c("-s", "--species"), help="choose the species <human>, <mouse> [default=%default]", default="human"),
#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
#make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
#make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-o", "--output"), help="additional tags for otuput [default=%default]", default="out")
#make_option(c("-f", "--fill_by"), help="choose the color you want to fill by [default=NA]", type='character', default=NA)
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)



##------------
## LIBRARIES
##------------ 

cat("Loading libraries... ")


suppressPackageStartupMessages(library("GO.db"))
if (opt$species == "human") {suppressPackageStartupMessages(library("org.Hs.eg.db"))}
if (opt$species == "mouse") {suppressPackageStartupMessages(library("org.Mm.eg.db"))}
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("plyr"))

cat("DONE\n\n")

############################
# BEGIN
############################

U = read.table(opt$universe, h=F, col.names='hs')
G = read.table(opt$genes, h=T, col.names='hs')


#orth = read.table('~/Documents/db/human-mouse/orthologs/merged_orthologs_HumanMouse_1_1_genes.tsv', h=F, col.names=c('hs','mm'))
#genes = read.table('~/Documents/human-mouse/tissues/maxdiff_tissuespecific.tsv', h=F, col.names=c('hs','mm','sample'))

# I want to create a list of parameters to perform GO enrichment on different gene sets

# take the entrez gene ids for all the orthologous genes which will be my universe (the same for all the sets)
if (opt$species == "human") {
universe = unlist(mget(U$hs, org.Hs.egENSEMBL2EG, ifnotfound=NA))}

if (opt$species == "mouse") {
universe = unlist(mget(U$hs, org.Mm.egENSEMBL2EG, ifnotfound=NA))}


sprintf("%s background genes; %s with a corresponding entrez id", nrow(U), length(unique(universe)))
# how many genes am I able to map?
# First thing notice that also ensembl gene ids longer than 15 characters are included
# if I remove these genes I end up with:
# length(unique(as.character(universe[which(nchar(names(universe)) == 15)]))) ----> 15593


createParams = function(x, species="human") {
	if (species == "human") {
	ann = "org.Hs.eg.db"
	geneset = unlist(mget(x, org.Hs.egENSEMBL2EG, ifnotfound=NA))}
	if (species == "mouse") {
	ann = "org.Mm.eg.db"
	geneset = unlist(mget(x, org.Mm.egENSEMBL2EG, ifnotfound=NA))}
	sprintf("%s foreground genes; %s with a corresponding entrez id", length(x), length(unique(geneset)))
	pv = 1-(1-0.05)**(1/length(x))
	params = new("GOHyperGParams",
		geneIds = geneset,
		universeGeneIds = universe,
		annotation = ann,
		ontology = opt$categ,
		pvalueCutoff = pv,
		conditional = TRUE,
		testDirection='over')
	return(params)}

res = hyperGTest(createParams(G$hs, opt$species))
write.table(summary(res), file=sprintf("%s.%s.tsv", opt$output, opt$categ), quote=F, sep="\t", row.names=F)
htmlReport(res, file=sprintf("%s.%s.html", opt$output, opt$categ))

q(save='no')


