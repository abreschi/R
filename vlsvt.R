#!/usr/bin/env Rscript

options(stringsAsFactors=F)

##------------
## LIBRARIES
##------------ 

suppressPackageStartupMessages(library("optparse"))

##################
# OPTION PARSING
##################

option_list <- list(

make_option(c("-i", "--input_matrix"), default="stdin", 
	help="the matrix you want to analyze. \"stdin\" for reading from standard input [default=%default]
		Columns are:
		col1 = transcript ID; col2 = gene ID; col3 = RPKM sample 1; col4 = RPKM sample 2; ...
		row = transcripts
	"),

make_option(c("-I", "--min.iso2.RPKM"), default=0.01,
	help="Minimum expression for the second major transcript. [default=%default]
		Particularly relevant to filter out genes expressing almost only one transcript. 
		Otherwise post-computation filtering is also possible.
	"),

make_option(c("-m", "--metadata"),  
	help="Matrix with metadata associated to each samples in the input matrix. If NULL vlsvt is computed across all samples [default=%default]"),

make_option(c("--merge_mdata_on"), default="labExpId",
	help="metadata column with headers of the input matrix [default=%default]"),

make_option(c("-f", "--factors"),
	help="column of the metadata matrix with the factor defining the groups"),

make_option(c("-o", "--output"), default="stdout",
	help="output file name. Can be stdout [default=%default]"),

make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
	help="Verbose output [default=%default]")

)

parser <- OptionParser(
	usage = "%prog [options] file", 
	option_list=option_list, 
	description="Compute VlsVt estimator across from RPKM values across an entire data.frame or between groups."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}

##   varGpGE: a data.frame with gene name, var.bwGp(i.e. variance explained by group classification), var.bwGp(i.e. gene expression contribution in variance explained by group classification), iso2.RPKM(i.e. the average RPKM expression of the second major isoform, for potential post-filtering).
##   groups, rescale: parameters used for the computation.


# ===============
# Functions
# ===============

# >>>>>>>>>>>>>>> VlsVt <<<<<<<<<<<<<<<<<<<<<

VlsVt <- function(data,rescale=TRUE){
  tol <- 10^-14
  k   <- dim(data)[2] # number of individuals
  n   <- dim(data)[1] # number of transcripts
  unC <- (k-1) / k    # to uncorrect
  if(rescale)
    data = sqrt(data)
  sigma <- unC * cov(t(data)) # corrected covariance
  ## Vt: Total variance
  Vt = sum(diag(sigma))
  ## Vls: variance in the model of constant splicing
  conjecture    <- svd(data,1,0)$u[,1]   
  indexP   <- conjecture > tol
  indexN   <- conjecture < -tol    
  nPos     <- sum(indexP)
  nNeg     <- sum(indexN)    
  if ((nPos > 0) & (nNeg>0)) { # Perron-Frobenius fails?
    require(NMFN)
    matrixFactor <- nnmf(data,1,method = "nnmf_als")
    bestDirection <- matrixFactor$W /  sqrt((t(matrixFactor$W) %*% matrixFactor$W)[1,1])
    Vls = ( t(bestDirection) %*% sigma ) %*% bestDirection
  }
  else {
    if (nNeg > 0) conjecture <- -1*conjecture 
    Vls = ( t(conjecture) %*% sigma ) %*% conjecture
  }
  ## Vls2: variance in the model of constant gene expression
  vecModel = vec = as.numeric(unlist(data))       # vectorize
  ## directament com OLS a veure si hi ha sort
  bigSum <- sum(vec)  
  for (i in 1:k) {
    littleSum                   <- k*sum(vec[(1+(i-1)*n):(i*n)])  
    vecModel[(1+(i-1)*n):(i*n)] <- vec[(1+(i-1)*n):(i*n)] + (bigSum - littleSum)/(n*k)
  }  
  nCompNeg    <- sum(vecModel < 0)
  if (nCompNeg == 0) data2g <- matrix(vecModel,nrow=n,ncol=k) 
  if (nCompNeg > 0) {
    Dg    <- matrix(0,nrow=(n*k),ncol=(n*k))
    Amatg <- matrix(0,nrow=(n*k+k-1),ncol=n*k)
    diag(Dg) = 1
    Amatg[1:(k-1),1:n] <- 1
    for (i in 1:(n*k)) {
      Amatg[(k-1+i),i] <- 1
      }
    for (i in 2:k) {
      Amatg[(i-1),(((i-1)*n)+1):(i*n)] <- -1
    }
    Amatg <- t(Amatg)
    dg    <- vec
    suppressPackageStartupMessages(library("quadprog"))
    sQP  <- solve.QP(Dg,dg,Amatg,meq=(k-1))
    Y2g         <- sQP$solution
    data2g       <- matrix(Y2g,nrow=n,ncol=k)
  }
  sigma2g    <- unC * cov(t(data2g)) ## Uncorrected variance
  trSigma2g  <- sum(diag(sigma2g))
  Vls2 = trSigma2g
  return(list(vls.vt=Vls/Vt,vls2.vt=Vls2/Vt))
}


# >>>>>>>>>>>>>>> VlsVt.bootstrap <<<<<<<<<<<<<<<<<<<<<

VlsVt.bootstrap <- function(data,rescale=TRUE,nb.bs=100,prop.bs=.5){
  if(!is.null(nb.bs)){
    require(plyr)
    vls.bs = ldply(1:nb.bs,function(ii){
      dat = data[,sample(1:ncol(data),ncol(data)*prop.bs,TRUE)]
      vlsvt.g = VlsVt(dat,rescale=rescale)
      data.frame(vls.vt = vlsvt.g$vls.vt,
                 vls2.vt = vlsvt.g$vls2.vt)
    })
    return(list(vls.vt=median(vls.bs$vls.vt),vls2.vt=median(vls.bs$vls2.vt)))
  } else {
    return(VlsVt(data,rescale=rescale))
  }
}


# >>>>>>>>>>>>>>> varGpGE <<<<<<<<<<<<<<<<<<<<<

varGpGE <- function(data,groups,rescale=TRUE){
    ##  'data' : the transcript expression data.
    ##  'group' : groups list: list of the column names of each group
  tol <- 10^-14
  data = data[,unlist(groups)]
  k   <- dim(data)[2] # number of individuals
  n   <- dim(data)[1] # number of transcripts
  ##
  if(rescale)
      data = sqrt(data)
  sigma <- cov(t(data)) # corrected covariance
  sigma = sigma * (k-1)
  sigma.bw = sigma
  for(gp in groups){
    sigma.gp = cov(t(data[,gp]))
    sigma.gp = sigma.gp * (length(gp)-1)
    sigma.bw = sigma.bw - sigma.gp
  }
  ## Vt
  ss.t = sum(diag(sigma))
  ss.t.bw = sum(diag(sigma.bw))
  var.tot.bw = ss.t.bw / ss.t
  ## Vls
  d.svd = svd(data,1,1)
  conjecture    <- d.svd$u[,1]
  indexP   <- conjecture > tol
  indexN   <- conjecture < -tol    
  nPos     <- sum(indexP)
  nNeg     <- sum(indexN)    
  if ((nPos > 0) & (nNeg>0)) { # Perron-Frobenius fails?
      ss.ge.bw = NA
  }  else {
    v = d.svd$v[,1]
    lambda = d.svd$d[1] * sum(conjecture) * as.numeric(v)
    dir = conjecture / sum(conjecture)
    data.p = lambda %*% t(dir)
    data.p = t(data.p)
    colnames(data.p) = colnames(data)
    sigma.p <- cov(t(data.p)) # corrected covariance
    sigma.p = sigma.p * (k-1)
    sigma.p.bw = sigma.p
    for(gp in groups){
      sigma.p.gp = cov(t(data.p[,gp]))
      sigma.p.gp = sigma.p.gp * (length(gp)-1)
      sigma.p.bw = sigma.p.bw - sigma.p.gp
    }
    ss.ge.bw = sum(diag(sigma.p.bw))
  }
  var.tot.bw.ge = ss.ge.bw / ss.t.bw
  return(list(var.tot.bw=var.tot.bw,
              var.tot.bw.ge=var.tot.bw.ge))
}


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

##   vls.vt: a data.frame with gene name, vls.vt(i.e. gene expression contribution estimate), vls2.vt(i.e. splicing contribution extimate), iso2.RPKM(i.e. the average RPKM expression of the second major isoform, for potential post-filtering).

# >>>>>>>>>>> VlsVt.comp <<<<<<<<<<<<<<<<

VlsVt.comp <- function(data,min.iso2.RPKM=opt$min.iso2.RPKM,rescale=TRUE,nb.bs=100,prop.bs=.5,verbose=opt$verbose){
  suppressPackageStartupMessages(library("plyr"))
  genes = unique(as.character(data[,2]))
  res = ldply(genes,function(gene){
    if(verbose) cat(gene,"...\n")
    if(length((gene.ids = which(data[,2]==gene)))>1){ ## At least two isoforms
      dat.g = data[gene.ids,-(1:2)]
      iso.RPKM.mean = rowMeans(dat.g,na.rm=TRUE)
	  iso1.RPKM = sort(iso.RPKM.mean,decreasing=TRUE)[1]
      iso2.RPKM = sort(iso.RPKM.mean,decreasing=TRUE)[2]
      vls.vt = vls2.vt = NA
      if(iso2.RPKM >= min.iso2.RPKM){ ## Second major isoform expressed
        vls.vt.o = VlsVt.bootstrap(dat.g,rescale=rescale,nb.bs=nb.bs,prop.bs=prop.bs)
        vls.vt = vls.vt.o$vls.vt
        vls2.vt = vls.vt.o$vls2.vt
      }
      return(data.frame(gene=gene,
                        vls.vt=vls.vt,
                        vls2.vt=vls2.vt,
                        iso2.RPKM=iso2.RPKM))
    }
  })
  return(res)
#  return(list(vls.vt=res,rescale=rescale,nb.bs=nb.bs,prop.bs=prop.bs))
}

#VlsVt.comp.multicore <- function(data,min.iso2.RPKM=.01,rescale=TRUE,nb.bs=100,prop.bs=.5,verbose=FALSE,nb.cores=4){
#  require(plyr)
#  require(parallel)
#  genes = unique(as.character(data[,2]))
#  res = mclapply(genes,function(gene){
#    tryCatch({
#      if(verbose) cat(gene,"...\n")
#      if(length((gene.ids = which(data[,2]==gene)))>1){ ## At least two isoforms
#        dat.g = data[gene.ids,-(1:2)]
#        iso.RPKM.mean = rowMeans(dat.g,na.rm=TRUE)
#        iso2.RPKM = sort(iso.RPKM.mean,decreasing=TRUE)[2]
#        vls.vt = vls2.vt = NA
#        if(iso2.RPKM >= min.iso2.RPKM){ ## Second major isoform expressed
#          vls.vt.o = VlsVt.bootstrap(dat.g,rescale=rescale,nb.bs=nb.bs,prop.bs=prop.bs)
#          vls.vt = vls.vt.o$vls.vt
#          vls2.vt = vls.vt.o$vls2.vt
#        }
#        return(data.frame(gene=gene,
#                          vls.vt=vls.vt,
#                          vls2.vt=vls2.vt,
#                          iso2.RPKM=iso2.RPKM))
#      } else {
#        return(data.frame(gene=gene,
#                          vls.vt=NA,
#                          vls2.vt=NA,
#                          iso2.RPKM=NA))
#      }}, error=function(e){
#        return(data.frame(gene=gene,
#                          vls.vt=NA,
#                          vls2.vt=NA,
#                          iso2.RPKM=NA))
#      })
#  },mc.cores=nb.cores)
#  res = ldply(res,identity)
#  return(list(vls.vt=res,rescale=rescale,nb.bs=nb.bs,prop.bs=prop.bs))
#}


# >>>>>>>>>>> varGpGE.comp <<<<<<<<<<<<<<<<

varGpGE.comp <- function(data,groups, min.iso2.RPKM=opt$min.iso2.RPKM, rescale=TRUE, verbose=opt$verbose){
  suppressPackageStartupMessages(library("plyr"))
  genes = unique(as.character(data[,2]))
  res = ldply(genes,function(gene){
    if(verbose) cat(gene,"...\n")
    if(length((gene.ids = which(data[,2]==gene)))>1){ ## At least two isoforms
      dat.g = data[gene.ids,-(1:2)]
      iso.RPKM.mean = rowMeans(dat.g,na.rm=TRUE)
      iso2.RPKM = sort(iso.RPKM.mean,decreasing=TRUE)[2]
      var.bwGp = var.bwGp.ge = NA
      if(iso2.RPKM >= min.iso2.RPKM){ ## Second major isoform expressed
          varGpGE.o = varGpGE(dat.g,groups,rescale=rescale)
          var.bwGp = varGpGE.o$var.tot.bw
          var.bwGp.ge = varGpGE.o$var.tot.bw.ge
      }
      return(data.frame(gene=gene,
                        var.bwGp=var.bwGp,
                        var.bwGp.ge=var.bwGp.ge,
                        iso2.RPKM=iso2.RPKM))
    }
  })
  return(res)
#  return(list(varGpGE=res,groups=groups,rescale=rescale))
}



###############
# BEGIN
###############

# Read matrix
if (opt$input_matrix == "stdin") {inF=file("stdin")} else {inF=opt$input_matrix}
m = read.table(inF, h=T, sep="\t")

# Read the metadata
if (!is.null(opt$metadata)) {
	mdata = read.table(opt$metadata, h=T, sep="\t", quote=NULL)
	mdata[,opt$merge_mdata_on] = gsub(",", ".", mdata[,opt$merge_mdata_on])
	mdata_col = unique(c(opt$merge_mdata_on, strsplit(opt$factors, "[+*:]")[[1]]))
	mdata = unique(mdata[,mdata_col])
	# Intersect mdata id with header of matrix
#	print(colnames(m)[-(1:2)] %in% mdata[,opt$merge_mdata_on])
	if (!all(colnames(m)[-(1:2)] %in% mdata[,opt$merge_mdata_on])) {
		cat("\n\tERROR: Some column names of the matrix are missing in the metadata\n\n")
		q(save='no')
	}
	mdata = mdata[mdata[,opt$merge_mdata_on] %in% colnames(m)[-(1:2)],]
	if (opt$verbose) {
	        cat("Metadata sample:\n")
	        print(head(mdata))
	}
	
	# Make the group list
	groups = lapply(split(mdata, mdata[,opt$factors]), function(x) x[,opt$merge_mdata_on])
#	print(groups)

	# Compute variance
	res = varGpGE.comp(m, groups)
} 

if (is.null(opt$metadata)) {
	if (opt$verbose) {
		cat("Across all samples...\n")
	}
	# Compute variance
	res = VlsVt.comp(m)
}


# ----------- Output ----------------

outF = ifelse(opt$output == "stdout", "", opt$output)
write.table(res, file=outF, sep="\t", quote=FALSE, row.names=FALSE)

q(save='no')
