library(RTCGA)
library(RTCGA.clinical)
infoTCGA()

dim(BRCA.clinical)
names(BRCA.clinical)
View(BRCA.clinical)

survivalTCGA(BRCA.clinical, OV.clinical, extract.cols = "admin.disease_code") -> BRCAOV.survInfo
kmTCGA(BRCAOV.survInfo, explanatory.names = "admin.disease_code",  pval = TRUE)

library(edgeR)
setwd("C:/Users/Srithegreat/Downloads/deseq_edger_test/")
rawdata <- read.delim("table_expression.txt", check.names=FALSE, stringsAsFactors=FALSE)
head(rawdata)
y <- DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])
library(org.Hs.eg.db)
idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
y <- y[idfound,]
egREFSEQ <- toTable(org.Hs.egREFSEQ)
m <- match(y$genes$RefSeqID, egREFSEQ$accession)
y$genes$EntrezGene <- egREFSEQ$gene_id[m]
y$genes$EntrezGene
egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol <- egSYMBOL$symbol[m]
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Symbol)
y <- y[!d,]
nrow(y)
y$samples$lib.size <- colSums(y$counts)
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL
y <- calcNormFactors(y)
plotMDS(y)
Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
data.frame(Sample=colnames(y),Patient,Tissue)
design <- model.matrix(~Patient+Tissue)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
#5% FDR
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
