#path to current working directory
setwd("/home/iob/Downloads/results_amit/results_amit/result_files")
inp<-read.csv(file="diffexpr-results_TH9_RA_TH9.csv")
head(inp)
library(ggplot2)
#ggrepel helps create non-overlapping labels
library(ggrepel)

#read the file with expression fold changes
genes<-read.csv("~/Desktop/diffexpr-results_TH9_TH0_new.csv")
head(genes)
#filter significant genes
genes$Significant <- ifelse(genes$padj< 0.05 & abs(genes$log2FoldChange)>1 , "Significant", "Not Sig")
#replace zero with least machine number precision
genes$pvalue<-replace(genes$pvalue,genes$pvalue==0,.Machine$double.eps)
ggplot(genes, aes(x = log2FoldChange, y = -log10(pvalue), text=Gene.Symbol)) +
  geom_point(aes(color = Significant), size=0.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")  +geom_text_repel(
    data = subset(genes,Highlight=="Yes"),
    aes(label = Gene.Symbol),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
#save in desired format
ggsave("volcano-1.png", width = 12, height = 8, dpi = 300)

#use plot.new to visualize using identify or setxy
plot.new()
plot(np)
#interactively identify points on graph
identify(x=genes$log2FoldChange, y = -log10(genes$pvalue), labels = genes$Gene.Symbol,pos=TRUE)
?identify

