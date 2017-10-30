
library(tximport)
library(gplots)
library(DESeq2)
library(readr)
setwd("C:/Users/Srithegreat/Downloads/results_amit/result_files/th9ra_th9/")
#read data into hash map entrezid,ensemblid
h<-read.table("C:/Users/Srithegreat/Downloads/entrezgene_ensemblmouse.txt",sep="\t")
mapping <- setNames(h[,1], h[,2])




samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samplesres
files <- file.path(dir, "rsem", samples$run, paste0(samples$run, ".genes.results"))
files<-list.files(".",pattern = "*.genes.results")#current directory files
files
names(files) <- c("TH0", "TH9_RA","TH9_RA","TH9_RA","TGFB","TGFB","TGFB","TH0","TH0","TH9","TH9","TH9","TH0_RA","TH0_RA","TH0_RA")
res1
head(txi.rsem$counts)
#import rsem using tximport
txi.rsem<-tximport(files,type="rsem")


#define conditions

sampleTable <- data.frame(condition = factor(c("TH0", "TH9_RA","TH9_RA","TH9_RA","TGFB","TGFB","TGFB","TH0","TH0","TH9","TH9","TH9","TH0_RA","TH0_RA","TH0_RA")))
sampleTable <- data.frame(condition = factor(c("TH9_RA", "TH9_RA","TH9_RA","TH9","TH9","TH9")))
rownames(sampleTable)<-c("1_1_1_2_rsem.genes.results","10_1_10_2_rsem.genes.results","11_1_11_2_rsem.genes.results","12_1_12_2_rsem.genes.results","13_1_13_2_rsem.genes.results","14_1_14_2_rsem.genes.results","15_1_15_2_rsem.genes.results","2_1_2_2_rsem.genes.results","3_1_3_2_rsem.genes.results","4_1_4_2_rsem.genes.results","5_1_5_2_rsem.genes.results","6_1_6_2_rsem.genes.results","7_1_7_2_rsem.genes.results","8_1_8_2_rsem.genes.results","9_1_9_2_rsem.genes.results")
rownames(sampleTable)<-c("1_1_1_2_rsem.genes.results","10_1_10_2_rsem.genes.results","11_1_11_2_rsem.genes.results","12_1_12_2_rsem.genes.results","2_1_2_2_rsem.genes.results","3_1_3_2_rsem.genes.results","4_1_4_2_rsem.genes.results","5_1_5_2_rsem.genes.results","6_1_6_2_rsem.genes.results")
colnames(txi.rsem$counts)
#condition<-factor(rep(c("TH0", "TH9","TH0_RA","TH9_RA","TGFB"), each = 3))
rownames(sampleTable) <- colnames(txi.rsem$counts)
colnames(txi.rsem$counts)
#run this next
colnames(txi.rsem$counts)<-rownames(sampleTable)
#replace any zero length transcripts with 1
txi.rsem$length<-replace(txi.rsem$length,txi.rsem$length==0,1)
rownames(sampleTable)
#create object for deseq2
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, design =~condition )
files
# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

#log transformation
rld <- rlogTransformation(dds)

#heatmap
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(sampleTable$condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(ggplot2)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[sampleTable$condition], RowSideColors=mycols[sampleTable$condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()
colData(rld)
#PCA custom
pcadata<-DESeq2::plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))
pcadata
ggplot(pcadata, aes(PC1, PC2, color=condition,shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

DESeq2::plotPCA(rld, intgroup="condition")
dev.new()


#specify the reference
dds$condition <- relevel(dds$condition, ref = "TH0")
resultsNames(dds)#check condition names here

# Get differential expression results
res <- results(dds,contrast=c("condition","TH9_RA","TH9"), alpha = 0.05)#Condition 1
res1 <- results(dds,contrast=c("condition","TH9","TH0"),alpha = 0.05)#Condition 2
res2 <- results(dds,contrast=c("condition","TH9_RA","TH0"),alpha = 0.05)#Condition 3
res3 <- results(dds,contrast=c("condition","TH0_RA","TH0"),alpha = 0.05)#COndition 4
res4 <- results(dds,contrast=c("condition","TH0_RA","TH9"),alpha = 0.05)#COndition 5
res5 <- results(dds,contrast=c("condition","TH9_RA","TH9"),alpha = 0.05)#COndition 6
res6 <- results(dds,contrast=c("condition","TGFB","TH9"),alpha = 0.05)#COndition 7
res7 <- results(dds,contrast=c("condition","TGFB","TH0_RA"),alpha = 0.05)#COndition 8
res8<-results(dds,contrast=c("condition","TH9_RA","TH0_RA"),alpha = 0.05)#COndition 9
res9 <- results(dds,contrast=c("condition","TGFB","TH9_RA"),alpha = 0.05)#COndition 10
?results
sum(res$padj<0.05,na.rm = TRUE)
table(res$padj<0.05 & abs(res$log2FoldChange)>1)
s## Order by adjusted p-value
res1 <- res1[order(res$padj), ]
res1

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=rownames(res), cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
#png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(res1, lfcthresh=4, sigthresh=0.01, textcx=.8, xlim=c(-2.3, 2))
dev.off()

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")
#plotMA
DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
rld
colnames(rld)<-c("TH0_1", "TH9_RA_1","TH9_RA_2","TH9_RA_3","TH0_2","TH0_3","TH9_1","TH9_2","TH9_3")
#adjusted p-values counts
sum(res$padj < 0.1, na.rm=TRUE)
#Heatmap of interesting genes
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )
f<-assay(rld)[ topVarGenes, ]
assay(rld)[topVarGenes,]
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), margin=c(8,12), cexRow = 0.7)

pheatmap(assay(rld)[ topVarGenes, ], 
         fontsize = 6, scale="row", dendrogram="column", color = colorRampPalette(rev(brewer.pal(n = 9, name =
                                                                                                  "RdYlGn")))(100))
nv<-add.anns(f,mart=ensembl)
nv
row.names(f)<-gsub("\\..*","",row.names(f))
#merge count data and print differentials
write.csv(nv,"annotated.csv")
resSig<-subset(res2,padj<=0.05)
resSig <- resSig[order(resSig$padj),]
#resdata <- merge(as.data.frame(resSig), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
#names(resSig)[1] <- "Gene"
#head(resdata)
row.names(resSig)<-gsub("\\..*","",row.names(resSig))
#annotate using Ensembl biomart
out<-add.anns(resSig,mart=ensembl)
row.names(out)<-NULL#remove first column for brevity
write.csv(out, file="diffexpr-results_TH9_RA_TH0.csv")


#interectively identify the differentials on graph
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]


# plot one gene of interest the read counts
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")


#visualizaion of results
library("ReportingTools")
library("org.Mm.eg.db")
library("biomaRt")
xx <- as.list(org.Mm.egENSEMBL2EG)
tmp=gsub("\\..*","",row.names(res))
ids<-xx[tmp]
#write.table(as.data.frame(resSig),"TGFB_TH0.tsv", sep="\t",quote=F)
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl" )
desReport <- HTMLReport(shortName = 'RNAseq_DESEQ2',
              title = 'RNA-seq analysis of differential expression using DESeq2',
                        reportDirectory = "./reports")
publish("TGFB vs TH0",desReport)
publish(res,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
          expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)
finish(desReport)

#map entrez IDs
entrezid<-mapIds(org.Mm.eg.db,
                 keys=row.names(res),
                 column="ENTREZID",
                 keytype="ENSEMBL",
                 multiVals="first")



## annotate using biomaRt
## note this is slightly different from what Mike pointed you to, as we
## are calling the IDs 'Ensembl', and are using mgi_symbol instead of hgnc_symbol
add.anns <- function(df, mart, ...) {
  nm <- rownames(df)
  anns <- getBM( attributes = c("ensembl_gene_id", "mgi_symbol", "description"), filters = "ensembl_gene_id", values = nm, mart = mart)
  anns <- anns[match(nm, anns[, 1]), ]
  colnames(anns) <- c("Ensembl", "Gene Symbol", "Gene Description")
  df <- cbind(anns, df[, 2:ncol(df)])
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
## assuming here that you have already created the test.dse object as before
temp <- HTMLReport("index","whatever")
publish(res3, desReport, .modifyDF = list(add.anns, modifyReportDF,ensemblLinks), make.plots=FALSE,mart = ensembl, pvalueCutoff = 0.05, factor = colData(dds)$condition)
## note that you have to use a large pvalueCutoff, because test.dse has large p-values
finish(desReport)
browseURL("index.html")
library(hwriter)


#publish all in one document
desReport <- HTMLReport(shortName = 'RNAseq_DESEQ2',
                        title = 'RNA-seq analysis of differential expression using DESeq2',
                        reportDirectory = "./reports")
publish("TGFB vs TH0",desReport)
publish(res,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)

publish("TH9 vs TH0",desReport)
publish(res1,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)
publish("TH9_RA vs TH0",desReport)
publish(res2,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)
publish("TH0_RA vs TH0",desReport)
publish(res3,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)
publish("TH0_RA vs TH9",desReport)
publish(res4,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)

publish("TH9_RA vs TH9",desReport)
publish(res5,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)
publish("TGFB vs TH9",desReport)
publish(res6,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)

publish("TGFB vs TH0_RA",desReport)
publish(res7,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)

publish("TH9_RA vs TH0_RA",desReport)
publish(res8,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)
publish("TGFB vs TH9_RA",desReport)
publish(res9,desReport,pvalueCutoff=0.05,
        factor = colData(dds)$condition,
        expName="deseq2",reportDir="./reports", make.plots=FALSE,.modifyDF = list(add.anns, modifyReportDF,ensemblLinks), mart = ensembl)



finish(desReport)




#new custom volcano

library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels
library("dplyr")
volcanoData

mutateddf <- mutate(as.data.frame(res1), sig=ifelse(res1$padj<0.05, "padj<0.05", "Not Sig")) #Will have different colors depending on significance
p = ggplot(mutateddf, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))
p
p+geom_text_repel(data=data=head(mutateddf, 20), aes(label=gene))





