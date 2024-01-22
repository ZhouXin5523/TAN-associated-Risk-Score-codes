# Fig.5 A
library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)

pFilter=0.001                
expFile="TCGA.TPM.txt"       
riskFile="risk.TCGA.txt"     
allDrugs=c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")


rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]


group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)


riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

for(drug in allDrugs){
	
	senstivity=pRRopheticPredict(data, drug, selection=1)
	senstivity=senstivity[senstivity!="NaN"]
	#senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
	
	
	sameSample=intersect(row.names(riskRT), names(senstivity))
	risk=riskRT[sameSample, "risk",drop=F]
	senstivity=senstivity[sameSample]
	rt=cbind(risk, senstivity)
	

	rt$risk=factor(rt$risk, levels=c("low", "high"))
	type=levels(factor(rt[,"risk"]))
	comp=combn(type, 2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
	
	test=wilcox.test(senstivity~risk, data=rt)
	
	if(test$p.value<pFilter){
		
		boxplot=ggboxplot(rt, x="risk", y="senstivity", fill="risk",
					      xlab="Risk",
					      ylab=paste0(drug, " senstivity (IC50)"),
					      legend.title="Risk",
					      palette=c("#0066FF","#FF0000")
					     )+ 
			stat_compare_means(comparisons=my_comparisons)
		pdf(file=paste0("durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
		print(boxplot)
		dev.off()
	}
}


# Fig.5 B
library(ggpubr)              
tciaFile="TCIA.txt"         
riskFile="risk.TCGA.txt"     



ips=read.table(tciaFile, header=T, sep="\t", check.names=F, row.names=1)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	

sameSample=intersect(row.names(ips), row.names(risk))
ips=ips[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(ips, risk)
	

data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


for(i in colnames(data)[1:(ncol(data)-1)]){
	rt=data[,c(i, "risk")]
	colnames(rt)=c("IPS", "Risk")
	gg1=ggviolin(rt, x="Risk", y="IPS", fill = "Risk", 
	         xlab="", ylab=i,
	         legend.title="Risk",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	         #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	
	pdf(file=paste0(i, ".pdf"), width=6, height=5)
	print(gg1)
	dev.off()
}



# Fig.5 C
library(limma)
library(plyr)
library(ggplot2)
library(ggpubr)

tideFile="TIDE.csv"           
riskFile="risk.TCGA.txt"     
  


tide=read.csv(tideFile, header=T, sep=",", check.names=F, row.names=1)
tide$Responder=ifelse(tide$Responder=="True", "Responder", "Non-responder")


group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0, ,drop=F]
row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=cbind(id=row.names(tide), tide)


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$risk=factor(risk$risk, levels=c("low", "high"))


sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "riskScore", drop=F]
data=cbind(tide, risk)


type=levels(factor(data[,"Responder"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


boxplot=ggboxplot(data, x="Responder", y="riskScore", fill="Responder",
		          xlab="",
		          ylab="Risk score",
		          legend.title="Responder"
		          #palette=bioCol
		          )+ 
	stat_compare_means(comparisons=my_comparisons)


pdf(file="boxplot.pdf", width=5.5, height=4.5)
print(boxplot)
dev.off()