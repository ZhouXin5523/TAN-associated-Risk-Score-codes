# Fig.6 A
library(maftools)       



risk=read.table("risk.TCGA.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)


geneNum=20    
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]


ann_colors=list()
col=c("blue", "red")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col


pdf(file="low risk matool.pdf", width=6, height=6)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()


pdf(file="high risk matool.pdf", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()



# Fig. 6B
library(limma)
library(ggpubr)
library(reshape2)



scoreCor=function(riskFile=null, scoreFile=null, project=null){
 
	data=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
	
	data=t(data)
	
	risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)
	
	sameSample=intersect(row.names(data),row.names(risk))
	data=data[sameSample,,drop=F]
	risk=risk[sameSample,,drop=F]
	rt=cbind(data,risk[,c("riskScore","risk")])
	rt=rt[,-(ncol(rt)-1)]
	
	
	immCell=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
	          "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
	          "Tfh","Th1_cells","Th2_cells","TIL","Treg")
	rt1=rt[,c(immCell,"risk")]
	data=melt(rt1,id.vars=c("risk"))
	colnames(data)=c("Risk","Type","Score")
	data$Risk=factor(data$Risk, levels=c("low","high"))
	p=ggboxplot(data, x="Type", y="Score", color = "Risk",
	     ylab="Score",add = "none",xlab="",palette = c("blue","red") )
	p=p+rotate_x_text(50)
	p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	pdf(file=paste0(project,".immCell.pdf"), width=7, height=6)
	print(p)
	dev.off()
	
	immFunction=c("APC_co_inhibition","APC_co_stimulation","CCR",
	          "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
	          "MHC_class_I","Parainflammation","T_cell_co-inhibition",
	          "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
	rt1=rt[,c(immFunction,"risk")]
	data=melt(rt1,id.vars=c("risk"))
	colnames(data)=c("Risk","Type","Score")
	data$Risk=factor(data$Risk, levels=c("low","high"))
	p=ggboxplot(data, x="Type", y="Score", color = "Risk",
	     ylab="Score",add = "none",xlab="",palette = c("blue","red") )
	p=p+rotate_x_text(50)
	p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	pdf(file=paste0(project,".immFunction.pdf"), width=7, height=6)
	print(p)
	dev.off()
}

scoreCor(riskFile="risk.TCGA.txt", scoreFile="immScore.all.txt", project="TCGA")
scoreCor(riskFile="risk.GEO.txt", scoreFile="immScore.all.txt", project="GEO")
scoreCor(riskFile="risk all.txt", scoreFile="immScore.all.txt", project="All")


#Fig.6C
library(limma)
library(ggpubr)
riskFile="risk.TCGA.txt"                          
subtypeFile="Subtype_Immune_Model_Based.txt"      


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


subtype=read.table(subtypeFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(subtype)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(subtype))

sameSample=intersect(row.names(subtype), row.names(risk))
subtype=subtype[sameSample,]
subtype=gsub(".+Immune |\\)","",subtype)
risk=risk[sameSample,]
data=cbind(as.data.frame(risk), subtype)
typeTab=table(data$subtype)
typeName=names(typeTab[typeTab>3])
data=data[which(data[,"subtype"] %in% typeName),]


group=levels(factor(data$subtype))
data$subtype=factor(data$subtype, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggboxplot(data, x="subtype", y="riskScore", color="subtype",
		          xlab="",
		          ylab="Risk score",
		          legend.title="Immune subtype",
		          add = "jitter")+ 
	    stat_compare_means(comparisons = my_comparisons)

pdf(file="immSubtype TCGA.pdf", width=5.5, height=5)
print(boxplot)
dev.off()

#Fig. 6D-E
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
riskFile="risk.TCGA.txt"       
RNAssFile="StemnessScores_RNAexp_20170127.2.tsv"      
DNAssFile="StemnessScores_DNAmeth_20170210.tsv"      



risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)


RNAss=read.table(RNAssFile, header=T, sep="\t",check.names=F, row.names=1)
RNAss=t(RNAss[1,,drop=F])
rownames(RNAss)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(RNAss))
RNAss=avereps(RNAss)


DNAss=read.table(DNAssFile, header=T, sep="\t", check.names=F, row.names=1)
DNAss=t(DNAss[1,,drop=F])
rownames(DNAss)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(DNAss))
DNAss=avereps(DNAss)


sameSample=Reduce(intersect,list(row.names(risk),row.names(DNAss),row.names(RNAss)))
risk=risk[sameSample,"riskScore",drop=F]
RNAss=RNAss[sameSample,,drop=F]
DNAss=DNAss[sameSample,,drop=F]
data=cbind(RNAss, DNAss, risk)

xlab="riskScore"
ylab="RNAss"
outFile="RNAss TCGA.cor.pdf"
x=as.numeric(data[,xlab])
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
		  xlab("Risk score") + ylab(ylab)+
		  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))

pdf(file=outFile, width=5.2, height=5)
print(p2)
dev.off()


xlab="riskScore"
ylab="DNAss"
outFile="DNAss TCGA.cor.pdf"
x=as.numeric(data[,xlab])
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
		  xlab("Risk score") + ylab(ylab)+
		  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))

pdf(file=outFile, width=5.2, height=5)
print(p2)
dev.off()



#Fig.6 F-H
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="TCGA all and GEO merge.txt"      
riskFile="risk all.txt"       
geneFile="checkpoint gene.txt"       

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	

gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
data=avereps(data)
	
risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
rt1=rt1[,c(sameGene,"risk")]


sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
	if(sd(rt1[,i])<0.001){next}
	wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"risk"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene=c(sigGene, "risk")
rt1=rt1[,sigGene]

rt1=melt(rt1,id.vars=c("risk"))
colnames(rt1)=c("risk","Gene","Expression")
	
group=levels(factor(rt1$risk))
rt1$risk=factor(rt1$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
	
boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="risk",
				  xlab="",
				  ylab="Gene expression",
				  legend.title="Risk",
				  width=0.8,
				  palette = c("#0066FF", "#FF0000") )+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=risk),
	method="wilcox.test",
	symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
	
pdf(file="checkpoint all.diff.pdf", width=8, height=5)
print(boxplot)
dev.off()


#Fig. 6I-J
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

library(ggridges)
library(ggplot2)

gene="riskScore"                               
gmtFile="c2.cp.kegg.v7.5.1.symbols.gmt"          

rt=read.table("gene risk.txt",sep="\t",header=T,check.names=F)
gmt=read.gmt(gmtFile)
	
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

dataL=data[,(data[gene,]<=median(data[gene,]))]    
dataH=data[,(data[gene,]>median(data[gene,]))]     
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH/meanL)
logFC=sort(logFC,decreasing=T)

kk=GSEA(logFC,TERM2GENE=gmt, nPerm=100,pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file=paste0("KEGG",".txt"),sep="\t",quote=F,row.names = F)
	

termNum=10
if(nrow(kkTab)>=termNum){
	   gsearidge=ridgeplot(kk,fill = "pvalue",showCategory = termNum,core_enrichment = TRUE,
                    label_format = 30,      
                    orderBy = "NES",
                    decreasing = F)
					pdf(file=paste0("KEGG ridge",".pdf"),width=8,height=11)
	print(gsearidge)
	dev.off()
         
	}

