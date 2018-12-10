## math 227A final presentation - sc-RNAseq pagoda2
library(pagoda2)

## define helper functions 
# set names euqal to the values
sn <- function(x) { names(x) <- x; return(x); }
# filter cells based on the gene/molecule dependency
t.filter.for.valid.cells <- function(countMatrix,min.cell.size=500, max.cell.size=5e4,p.level=min(1e-3,1/ncol(countMatrix)),alpha=0.1,do.par=T) {
  if(do.par) { par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);}
  hist(log10(colSums(countMatrix)),col='wheat',xlab='log10[ molecules ]',main='') 
  # some of the cells are very large .. those can skew the analysis of more subtle populations (too much bias) .. letting them in here though
  
  abline(v=log10(c(min.cell.size,max.cell.size)),lty=2,col=2)
  # look at the number of genes vs. molecule size depenency
  df <- data.frame(molecules=colSums(countMatrix),genes=colSums(countMatrix>0)); 
  df <- df[df$molecules>=min.cell.size,];
  df <- log10(df);
  df <- df[order(df$molecules,decreasing=F),]
  plot(df,col=adjustcolor(1,alpha=alpha),cex=0.5,ylab='log10[ gene counts]',xlab='log10[ molecule counts]')
  abline(v=log10(c(min.cell.size,max.cell.size)),lty=2,col=2)
  #abline(lm(genes ~ molecules, data=df),col=4)
  require(MASS)  
  m <- rlm(genes~molecules,data=df)
  suppressWarnings(pb <- data.frame(predict(m,interval='prediction',level = 1-p.level,type="response")))
  polygon(c(df$molecules,rev(df$molecules)),c(pb$lwr,rev(pb$upr)),col=adjustcolor(2,alpha=0.1),border = NA)
  outliers <- rownames(df)[df$genes > pb$upr | df$genes < pb$lwr];
  points(df[outliers,],col=2,cex=0.6)
  # set of filtered cells to move forward with  
  valid.cells <- colSums(countMatrix)>min.cell.size & colSums(countMatrix)<max.cell.size & !(colnames(countMatrix) %in% outliers)
  countMatrix[,valid.cells,drop=F]
}
# load 10x matrices from a named list of result folders
t.load.10x.data <- function(matrixPaths) {
  require(parallel)
  require(Matrix)
  mclapply(sn(names(matrixPaths)),function(nam) {
    matrixPath <- matrixPaths[nam];
    # read all count files (*_unique.counts) under a given path
    #cat("loading data from ",matrixPath, " ");
    x <- as(readMM(paste(matrixPath,'matrix.mtx',sep='/')),'dgCMatrix'); # convert to the required sparse matrix representation
    cat(".")
    gs <- read.delim(paste(matrixPath,'genes.tsv',sep='/'),header=F)
    rownames(x) <- gs[,2]
    cat(".")
    gs <- read.delim(paste(matrixPath,'barcodes.tsv',sep='/'),header=F)
    colnames(x) <- gs[,1]
    cat(".")
    colnames(x) <- paste(nam,colnames(x),sep='_');
    x
  },mc.cores=30)
}


## plot inaccuracy bar plot for hand-writing digit
cd <- t.load.10x.data(list(PBMC8K='/Users/yanwengong/Documents/fall_2018/math_227A/final_presentation/raw_gene_bc_matrices/GRCh38'))
str(cd) ## dim:33694(gene)x737280(cell)
## require at least 500 molecules per cells
counts <- t.filter.for.valid.cells(cd[[1]],min.cell.size=500)
## chcek number of molecules per genesm and omit low-expression genes
hist(log10(rowSums(counts)+1),main='Molecules per gene',xlab='molecules (log10)',col='wheat')
abline(v=1,lty=2,col=2)
counts <- counts[rowSums(counts)>=10,]
str(counts) ## now it's 15556x8786

## make gene names unique and build pagoda object
rownames(counts) <- make.unique(rownames(counts))
r <- Pagoda2$new(counts,log.scale=FALSE)

## adjust variance 
r$adjustVariance(plot=T,gam.k=10)

## dimension reduction with PCA
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine');
r$getKnnClusters(method=infomap.community,type='PCA')

## k-nearest neighbor graph space for clustering and visualization calculations
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=F)
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='clusters (tSNE)')

## look at different depth and expression pattern of one gene
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$depth,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='depth')
gene <-"LYZ"
r$plotEmbedding(type='PCA',embeddingType='tSNE',colors=r$counts[,gene],shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main=gene)

## multiple clustering methods can be used
## FOR THIS FUNCTION, DO THEY USE 100 PCs to PERFORM CLUSTERING???
r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
par(mfrow=c(1,2))
r$plotEmbedding(type='PCA',embeddingType='tSNE',groups=r$clusters$PCA$community,show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='infomap clusters (tSNE)')
r$plotEmbedding(type='PCA',embeddingType='tSNE',clusterType='multilevel',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='multlevel clusters (tSNE)')
