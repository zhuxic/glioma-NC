#引用包
library(limma)
library(ComplexHeatmap)
expFile="geneExp.txt"       #表达数据文件
cliFile="GCN5clinical.txt"      #临床数据文件
setwd("H:\\shengxin\\cgga\\gGCN5\\g\\clinikheatp")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

#删掉正常样品
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

#根据目标基因表达量对样品进行分组
Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
Type=factor(Type, levels=c("Low","High"))
data=cbind(as.data.frame(data), Type)
data=data[order(data[,gene]),] 

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>41,">41","<=41"))

#合并数据
samSample=intersect(row.names(data), row.names(cli))
data=data[samSample,"Type",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(data, cli)

#对临床性状进行循环，观察临床性状在高低表达组之间是否具有差异，得到显著性标记
sigVec=c(gene)
for(clinical in colnames(rt[,2:ncol(rt)])){
	data=rt[c("Type", clinical)]
	colnames(data)=c("Type", "clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	tableStat=table(data)
	stat=chisq.test(tableStat)
	pvalue=stat$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(clinical, Sig))
}
colnames(rt)=sigVec

#定义热图注释的颜色
#rt=rt[apply(rt,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
bioCol=c("#0066FF","#FF9900","#FF0000","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
colorList=list()
colorList[[gene]]=c("Low"="blue", "High"="red")
j=0
for(cli in colnames(rt[,2:ncol(rt)])){
	cliLength=length(levels(factor(rt[,cli])))
	cliCol=bioCol[(j+1):(j+cliLength)]
	j=j+cliLength
	names(cliCol)=levels(factor(rt[,cli]))
	cliCol["unknow"]="grey75"
	colorList[[cli]]=cliCol
}

#绘制热图
ha=HeatmapAnnotation(df=rt, col=colorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(rt))
Hm=Heatmap(zero_row_mat, top_annotation=ha)

#输出热图
pdf(file="mergeheatmap2.pdf", width=12, height=13)
draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

####################################################################################
#install.packages("pheatmap")

gene="KAT2A"                   #按哪个基因进行分组
logFCfilter=1                   #logFC过滤阈值
adjPfilter=0.05                   #矫正后p值阈值

library(limma)                    #引用包
setwd("H:\\shengxin\\cgga\\gGCN5\\g")               #设置工作目录
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
pdf(file="heatmap2.pdf",height=6,width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 10,
         fontsize_row=8,
         fontsize_col=10)
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

