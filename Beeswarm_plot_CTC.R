f<-read.table("C:/Users/Srithegreat/Desktop/CTC_file.txt",sep="\t",header=T)
flog<-log10(f+1)
boxplot(f,lwd = 2, ylab = "Generic EMT Score", las=2, ylim=c(-1,1))
stripchart(f, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, col = ifelse(x > 0,'blue','green'))

legend("topleft", fill=c("lightblue","yellow"), legend=c("Cad7","Ecad"))
?legend
plot(density(f$X0secCAD))
density(f$X0secCAD)
hist(breaks)
library('beeswarm')
par(mfrow=c(2,2))
par(mfrow=c(1,1))
beeswarm(f,col=rainbow(5), method="square", pch=16, labels=colnames(f))
beeswarm(f,col=rainbow(5), method="hex", pch=16, labels=colnames(f))
beeswarm(f,col=rainbow(5), method="swarm", pch=16, labels=colnames(f))

beeswarm(f,col = ifelse(f > 0,'blue','green'), method="square", pch=16, labels=colnames(f), las=2, ylim=c(-1,1), ylab = "Generic EMT Score")
abline(h=0, lty=3)
install.packages("corrgram")
library(corrgram)
f<-read.table("C:/Users/Srithegreat/Downloads/heatmap_ccle.txt",header=TRUE,sep="\t")
d<-cor(f[,3:20])
scatter(d)
d<-corrgram(f[,3:20],order=NULL,lower.panel=panel.shade,
            upper.panel=NULL, text.panel=panel.txt)
heatmap(d,na.rm = TRUE, Rowv = FALSE,Colv = FALSE,hclustfun = "none")
d
