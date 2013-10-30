##------------
## LIBRARIES
##------------ 
cat("Loading libraries... ")
suppressPackageStartupMessages(library("PFAM.db"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))
cat("DONE\n\n")

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-u", "--universe"), help="a list of human gene identifiers (ensEMBL ids), NO header"),
make_option(c("-G", "--genes"), help="a list of human gene identifiers for the foreground (ensEMBL ids), WITH header"),
make_option(c("-o", "--output"), help="additional tags for otuput [default=out]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
print(opt)

na2null = function(x) if(is.na(x)) {return(NULL)}else{return(x)}


############################
# BEGIN
############################

U = read.table(opt$universe, h=F, col.names='hs')
G = read.table(opt$genes, h=T, col.names='hs')


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
	params = new("PFAMHyperGParams",
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
