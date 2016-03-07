#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("optparse"))
options(stringsAsFactors=F)

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-u", "--universe"), help="a list of human gene identifiers (ensEMBL ids), NO header"),
make_option(c("-G", "--genes"), help="a list of human gene identifiers for the foreground (ensEMBL ids), WITH header"),
#make_option(c("-l", "--log"), action="store_true", default=FALSE, help="apply the log [default=FALSE]"),
#make_option(c("-p", "--pseudocount"), type="double", help=sprintf("specify a pseudocount for the log [default=%s]",pseudocount), default=pseudocount),
#make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-o", "--output"), help="additional tags for otuput [default=out]")
#make_option(c("-f", "--fill_by"), help="choose the color you want to fill by [default=NA]", type='character', default=NA)
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)

##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library("KEGG.db"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("plyr"))

na2null = function(x) if(is.na(x)) {return(NULL)}else{return(x)}


############################
# BEGIN
############################

U = read.table(opt$universe, h=F, col.names='hs')
G = read.table(opt$genes, h=T, col.names='hs')


#orth = read.table('~/Documents/db/human-mouse/orthologs/merged_orthologs_HumanMouse_1_1_genes.tsv', h=F, col.names=c('hs','mm'))
#genes = read.table('~/Documents/human-mouse/tissues/maxdiff_tissuespecific.tsv', h=F, col.names=c('hs','mm','sample'))

# I want to create a list of parameters to perform GO enrichment on different gene sets

# take the entrez gene ids for all the orthologous genes which will be my universe (the same for all the sets)
universe = unlist(mget(U$hs, org.Hs.egENSEMBL2EG, ifnotfound=NA))

sprintf("%s background genes; %s with a corresponding entrez id", nrow(U), length(unique(universe)))
# how many genes am I able to map?
# First thing notice that also ensembl gene ids longer than 15 characters are included
# if I remove these genes I end up with:
# length(unique(as.character(universe[which(nchar(names(universe)) == 15)]))) ----> 15593


createParams = function(x) {
	geneset = unlist(mget(x, org.Hs.egENSEMBL2EG, ifnotfound=NA))
	sprintf("%s foreground genes; %s with a corresponding entrez id", length(x), length(unique(geneset)))
	pv = 1-(1-0.05)**(1/length(x))
	params = new("KEGGHyperGParams",
		geneIds = geneset,
		universeGeneIds = universe,
		annotation = 'org.Hs.eg.db',
		pvalueCutoff = pv,
		testDirection='over')
	return(params)}

res = hyperGTest(createParams(G$hs))
write.table(summary(res), file=sprintf("%s.tsv", opt$output), quote=F, sep="\t", row.names=F)
htmlReport(res, file=sprintf("%s.html", opt$output))

q(save='no')


#listOfParamObjs = aggregate(G$hs, list(G$sample), createParams)$x

# run it (and hope it is ok! )
#resultList = lapply(listOfParamObjs, hyperGTest)

# give the name to each GO enrichment
#names(resultList) <- aggregate(genes$hs, list(genes$sample), length)[,1]

# put the results in a data.frame
#resultDf = rbind.fill(lapply(as.list(names(resultList)[-which(lapply(resultList, function(x) nrow(summary(x))) == 0)]), 
#	function(n)data.frame(summary(resultList[[n]]), n_genes=length(geneIds(resultList[[n]])), sample=n)))

# save the data.frame in a file
#write.table(resultDf, file='maxdiff_tissuespecific_GO.tsv', quote=F, row.names=F, col.names=T, sep='\t')

