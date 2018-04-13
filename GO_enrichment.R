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
make_option(c("-u", "--universe"), help="a list of gene identifiers (ensEMBL ids), NO header"),
make_option(c("-G", "--genes"), default="stdin",
	help="a list of gene identifiers for the foreground (ensEMBL ids), WITHOUT header [default=%default]"),
make_option(c("-c", "--categ"), help="choose the GO category < BP | MF | CC > [default=%default]", default="BP"),
make_option(c("-s", "--species"), help="choose the species < human | mouse | dmel > [default=%default]", default="human"),
make_option(c("-o", "--output"), help="additional tags for otuput [default=%default]", default="out"),
make_option(c("-d", "--output_dir"), default="./", help="directory for the output [default=%default]"),
make_option(c("--output_genes"), action="store_true", default=FALSE, help="Output the list of genes from enriched GO terms in a separate file [default=%default]"),
make_option(c("-v", "--verbose"), default=FALSE, action="store_true", help="verbose")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}



##------------
## LIBRARIES
##------------ 

if (opt$verbose) {cat("Loading libraries... ")}


suppressPackageStartupMessages(library("GO.db"))
if (opt$species == "human") {ann = "org.Hs.eg.db"; suppressPackageStartupMessages(library("org.Hs.eg.db"))}
if (opt$species == "mouse") {ann = "org.Mm.eg.db"; suppressPackageStartupMessages(library("org.Mm.eg.db"))}
if (opt$species == "dmel") {ann = "org.Dm.eg.db"; suppressPackageStartupMessages(library("org.Dm.eg.db"))}
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("plyr"))

if (opt$verbose) {cat("DONE\n\n"); sessionInfo()}

############################
# BEGIN
############################

U = read.table(opt$universe, h=F, col.names='hs')

U$hs = unique(U$hs)

#if (opt$genes == "stdin") {
#	G = read.table(file("stdin"), h=T, col.names='hs')	
#} else {
#	G = read.table(opt$genes, h=T, col.names='hs')
#}

if (opt$genes == "stdin") {
	G = read.table(file("stdin"), h=F, col.names='hs')	
} else {
	G = read.table(opt$genes, h=F, col.names='hs')
}


# I want to create a list of parameters to perform GO enrichment on different gene sets

# take the entrez gene ids for all the orthologous genes which will be my universe (the same for all the sets)
if (opt$species == "human") {
universe = unlist(mget(U$hs, org.Hs.egENSEMBL2EG, ifnotfound=NA))}

if (opt$species == "mouse") {
universe = unlist(mget(U$hs, org.Mm.egENSEMBL2EG, ifnotfound=NA))}

if (opt$species == "dmel") {
universe = unlist(mget(U$hs, org.Dm.egENSEMBL2EG, ifnotfound=NA))}


if (opt$verbose) {sprintf("%s background genes; %s with a corresponding entrez id", nrow(U), length(unique(universe)))}
# how many genes am I able to map?
# First thing notice that also ensembl gene ids longer than 15 characters are included
# if I remove these genes I end up with:
# length(unique(as.character(universe[which(nchar(names(universe)) == 15)]))) ----> 15593


createParams = function(x, species="human") {
	if (species == "human") {
		geneset = unlist(mget(x, org.Hs.egENSEMBL2EG, ifnotfound=NA))
	}
	if (species == "mouse") {
		geneset = unlist(mget(x, org.Mm.egENSEMBL2EG, ifnotfound=NA))
	}
	if (species == "dmel") {
		geneset = unlist(mget(x, org.Dm.egENSEMBL2EG, ifnotfound=NA))
	}
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

res = hyperGTest(createParams(unique(G$hs), opt$species))

if (opt$verbose) {
cat ("Finished Hypergeometric test\n")
}

# Reformat the output table
df = summary(res)
df$Pvalue = format(df$Pvalue, digits=1) 
df$OddsRatio <- round(df$OddsRatio, 2)
df$ExpCount <- round(df$ExpCount, 2)



# Get the genes for the enriched GO terms
if (opt$output_genes) {
	enrichGenes <- ldply(geneIdsByCategory(res, catids=sigCategories(res, pvalueCutoff(res))), data.frame)
	enrichGenes[[2]] <- mapIds(eval(parse(text=ann)), keys=enrichGenes[[2]], keytype="ENTREZID", column=c("ENSEMBL"))
	colnames(enrichGenes) <- c("GO", "gene_id")
	output = sprintf("%s/%s.%s.genes", opt$output_dir, opt$output, opt$categ)
	write.table(enrichGenes, file=output, quote=F, sep="\t", row.names=F)
}


# Print output
output = sprintf("%s/%s.%s", opt$output_dir, opt$output, opt$categ)
write.table(df, file=sprintf("%s.tsv", output), quote=F, sep="\t", row.names=F)
htmlReport(res, file=sprintf("%s.html", output))




q(save='no')


