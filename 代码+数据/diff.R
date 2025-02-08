

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

gene="KAT2A"                   #按哪个基因进行分组
logFCfilter=1                   #logFC过滤阈值
adjPfilter=0.05                   #矫正后p值阈值

library(limma)                    #引用包
setwd("H:\\shengxin\\cgga\\gGCN5\\g\\GCN5hp")               #设置工作目录
rt=read.table("normalize.txt",sep="\t",header=T,check.names=F)         #读取输入文件

#如果一个基因存在多行，取均值
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

#按照基因表达中位值对样品分组
low=rt[gene,]<=median(rt[gene,])
high=rt[gene,]>median(rt[gene,])
lowRT=rt[,low]
highRT=rt[,high]
conNum=ncol(lowRT)        #低表达组样品数目
treatNum=ncol(highRT)     #高表达组样品数目
rt=cbind(lowRT,highRT)

#差异分析
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="all.xls",sep="\t",quote=F)

#输出差异结果
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adjPfilter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="diff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut,file="diff.txt",sep="\t",quote=F,col.names=F)

#绘制差异基因热图
library(pheatmap)
geneNum=20
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=rt[hmGene,]
Type=c(rep("Low",conNum),rep("High",treatNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=8,width=12)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 10,
         fontsize_row=12,
         fontsize_col=14)
dev.off()

#输出相关性分析的数据
hmExpOut=rbind(id=colnames(hmExp),hmExp)
write.table(hmExpOut,file="corInput.txt",sep="\t",quote=F,col.names=F)


#火山图
pdf(file="vol.pdf",width=5,height=5)
yMax=80
xMax=40
plot(allDiff$logFC, -log10(allDiff$adj.P.Val), ylab="-log10(adj.P.Val)",xlab="logFC",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1)
diffSub=subset(allDiff, adj.P.Val<adjPfilter & logFC>logFCfilter)
points(diffSub$logFC, -log10(diffSub$adj.P.Val), pch=20, col="red",cex=1.2)
diffSub=subset(allDiff, adj.P.Val<adjPfilter & logFC<(-logFCfilter))
points(diffSub$logFC, -log10(diffSub$adj.P.Val), pch=20, col="green",cex=1.2)
abline(v=0,lty=2,lwd=3)
dev.off()

