library(edgeR)
library(tximport)
library(biomaRt)
library(ReportingTools)
library(hwriter)
setwd("C:/Users/Srithegreat/Downloads/results_amit/result_files/")
files<-list.files(".")#current directory files

txi.rsem <- tximport(files, type = "rsem")
names(files)<-c("TH0", "TH9_RA","TH9_RA","TH9_RA","TGFB","TGFB",
                "TGFB","TH0","TH0","TH9","TH9","TH9","TH0_RA","TH0_RA","TH0_RA")

cts <- txi.rsem$counts
#replace zero length
txi.rsem$length<-replace(txi.rsem$length,txi.rsem$length==0,1)
normMat <- txi.rsem$length
normMat <- normMat/exp(rowMeans(log(normMat)))
#estimate normalization factors
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
condition = factor(c("TH0", "TH9_RA","TH9_RA","TH9_RA","TGFB","TGFB","TGFB","TH0","TH0","TH9","TH9","TH9","TH0_RA","TH0_RA","TH0_RA"))
y <- DGEList(cts,group = condition)
y$offset <- t(t(log(normMat)) + o)
#crete model
design <- model.matrix(~0+condition)
design
#estimate dispersion

d2 = estimateGLMTrendedDisp(y, design)
d2 = estimateGLMTagwiseDisp(d2, design)

#fit a GLM
f = glmFit(d2, design)
#perform likelihood test
#THO_TGFB
warnings()
#create contrast matrices for all comparisons
de = glmLRT(f,contrast=c(1,-1,0,0,0))
de1 = glmLRT(f,contrast=c(0,-1,1,0,0))
de2 = glmLRT(f,contrast=c(0,-1,0,1,0))
de3 = glmLRT(f,contrast=c(0,-1,0,0,1))
de4 = glmLRT(f,contrast=c(0,0,1,-1,0))
de5 = glmLRT(f,contrast=c(0,0,0,-1,1))
de6 = glmLRT(f,contrast=c(1,0,0,-1,0))
de7 = glmLRT(f,contrast=c(1,0,-1,0,0))
de8 = glmLRT(f,contrast=c(0,0,-1,0,1))
de9 = glmLRT(f,contrast=c(1,0,0,0,-1))
#differential expression test
#et <- exactTest(y, pair=c("TH9","TH0"))
options(digits=3)
tt<-topTags(de9,p.value=0.05,n=nrow(y))
#summary(decideTestsDGE(de))

#remove accession versions
rownames(tt$table)<-gsub("\\..*","",row.names(tt$table))
#annotate ensembl
out<-add.anns(tt$table,mart=ensembl)
row.names(out)<-NULL#remove first column for brevity
write.csv(out, file="toptags_TGFB_TH9RA.csv",row.names=FALSE)


ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl" )
desReport <- HTMLReport(shortName = 'RNAseq_EDGER',
                        title = 'RNA-seq analysis of differential expression using EDGER',
                        reportDirectory = "./reports")
publish("TGFB vs TH0",desReport)
publish(tt$table,desReport,pvalueCutoff=0.05,factor =condition ,
        expName="edgeR",reportDir="./reports", .modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)
finish(desReport)


#inspect relationship between samples
plotMDS(y, labels=samples$shortname,col=c("darkgreen","blue")[factor(samples$condition)])

#inspect depth adjusted reads
nc = cpm(y, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
head(nc[rn,order(condition)],5)
#Summarize
summary(de <- decideTestsDGE(de, p=0.05))
#plot M vs A 5% FDR
deg = rn[tt$table$FDR < .05]
plotSmear(y, de.tags=deg)
abline(h=c(-1,1), col="blue")

#write results
write.csv(tt$table, file="toptags_edgeR_TGFB_TH0.csv", row.names=TRUE)

#Perform GO
go <- goana(tt, species="Mm")
topGO(go, sort="up")
keg <- kegga(tt, species="mm")
topKEGG(keg, sort="up")



#annotate using ensembl
add.anns <- function(df, mart, ...) {
  nm <- rownames(df)
  anns <- getBM( attributes = c("ensembl_gene_id", "mgi_symbol", "description"), filters = "ensembl_gene_id", values = nm, mart = mart)
  anns <- anns[match(nm, anns[, 1]), ]
  colnames(anns) <- c("Ensembl", "Gene Symbol", "Gene Description")
  df <- cbind(anns, df[, 1:ncol(df)])
  rownames(df) <- nm
  df
}

#Add links to Ensembl.org, because that's how we roll.
ensemblLinks <- function(df, ...){
  naind <- is.na(df$Ensembl)
  df$Ensembl <- hwrite(as.character(df$Ensembl), link = paste0("http://www.ensembl.org/Mus_musculus/Gene/Summary?g=",
                                                               as.character(df$Ensembl)), table = FALSE)
  df$Ensembl[naind] <- ""
  return(df)
}

library(org.Mm.eg.db)

#map IDS and add annotation
tt$table$symbol <- mapIds(org.Mm.eg.db,
                     keys=row.names(tt$table),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
tt$table$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(tt$table),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
tt$table$description <- mapIds(org.Mm.eg.db,
                          keys=row.names(tt$table),
                          column="GENENAME",
                          keytype="ENSEMBL",
                          multiVals="first")



#export the results
allfiles<-list.files(".",pattern="toptags.*csv")
desReport <- HTMLReport(shortName = 'RNAseq_EDGER',
                        title = 'RNA-seq analysis of differential expression using EDGER',
                        reportDirectory = "./EDGER_RESULTS")
for (i in allfiles){
  
  inp<-read.csv(i,header=TRUE)
  #get partial name for comparision  
  f<-gsub(".csv","",gsub("toptags_","",i));
  publish(f,desReport)
  publish(inp,desReport,expName="edgeR",reportDir="./EDGER_RESULTS", .modifyDF = modifyReportDF)
  
}
finish(desReport)