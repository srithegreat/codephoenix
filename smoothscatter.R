fname<-read.table("C:/Users/Srithegreat/Downloads/1ug vs 2ug.txt",header=T,sep = "\t",row.names=1)
f<-hist(log2(fname$X2ug.1ug))
plot(f$breaks,f$counts type="o")
warnings()
axis(1,at=f$breaks,labels=f$breaks)
f$breaks
f$counts
n<-read.table("C:/Users/Srithegreat/Desktop/forcomparison_prot_trans.txt",header=T,sep="\t")
par(mfrow=c(1,2))
smoothScatter(n$log10protein,n$protein.mRNA)
legend("topleft",legend=paste("R=",round(cor(n$log10protein,n$protein.mRNA, method="spearman", use="na.or.complete"),1)),bty="n" )

?cor
par(mfrow=c(2,3))
smoothScatter(n$log10protein,n$EstBladder, main="Bladder", xlab="log10(ProteinIntensity)", ylab="Estimated_protein")
legend("topleft",legend=paste("corr=",round(cor(n$log10protein,n$EstBladder, method="spearman"),2)),bty="n" )

smoothScatter(n$log10protein,n$EstColon, main="Colon", xlab="log10(ProteinIntensity)", ylab="Estimated_protein")
legend("topleft",legend=paste("corr=",round(cor(n$log10protein,n$EstColon, method="spearman"),2)),bty="n" )

smoothScatter(n$log10protein,n$EstIleum, main="Ileum", xlab="log10(ProteinIntensity)", ylab="Estimated_protein")
legend("topleft",legend=paste("corr=",round(cor(n$log10protein,n$EstIleum, method="spearman"),2)),bty="n" )
smoothScatter(n$log10protein,n$EstProstate, main="Prostate", xlab="log10(ProteinIntensity)", ylab="Estimated_protein")
legend("topleft",legend=paste("corr=",round(cor(n$log10protein,n$EstProstate, method="spearman"),2)),bty="n" )
smoothScatter(n$log10protein,n$EstTestis, main="Testis", xlab="log10(ProteinIntensity)", ylab="Estimated_protein")
legend("topleft",legend=paste("corr=",round(cor(n$log10protein,n$EstTestis, method="spearman"),2)),bty="n" )
smoothScatter(n$log10protein,n$EstRectum, main="Rectum", xlab="log10(ProteinIntensity)", ylab="Estimated_protein")
legend("topleft",legend=paste("corr=",round(cor(n$log10protein,n$EstRectum, method="spearman"),2)),bty="n" )
