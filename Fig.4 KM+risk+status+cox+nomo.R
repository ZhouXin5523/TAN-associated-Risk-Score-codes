#Fig.4 A-C

install.packages("pheatmap")

library(pheatmap)           

bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
		rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)   
		rt=rt[order(rt$riskScore),]                                            
		
		
		riskClass=rt[,"risk"]
		lowLength=length(riskClass[riskClass=="low"])
		highLength=length(riskClass[riskClass=="high"])
		line=rt[,"riskScore"]
		line[line>10]=10
		pdf(file=riskScoreFile,width = 10,height = 3.5)
		plot(line, type="p", pch=20,
		     xlab="Patients (increasing risk socre)", ylab="Risk score",
		     col=c(rep("green",lowLength),rep("red",highLength)) )
		abline(h=median(rt$riskScore),v=lowLength,lty=2)
		legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
		dev.off()
		
		
		color=as.vector(rt$fustat)
		color[color==1]="red"
		color[color==0]="green"
		pdf(file=survStatFile,width = 10,height = 3.5)
		plot(rt$futime, pch=19,
		     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
		     col=color)
		legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
		abline(v=lowLength,lty=2)
		dev.off()
		
		
		rt1=rt[c(3:(ncol(rt)-2))]
		rt1=log2(rt1+1)
		rt1=t(rt1)
		annotation=data.frame(type=rt[,ncol(rt)])
		rownames(annotation)=rownames(rt)
		pdf(file=heatmapFile,width = 10,height = 3.5)
		pheatmap(rt1, 
		         annotation=annotation, 
		         cluster_cols = FALSE,
		         fontsize_row=11,
		         show_colnames = F,
		         fontsize_col=3,
		         color = colorRampPalette(c("green", "black", "red"))(50) )
		dev.off()
}
bioRiskPlot(inputFile="risk.GEO.txt",riskScoreFile="geo.riskScore.pdf",survStatFile="geo.survStat.pdf",heatmapFile="geo.heatmap.pdf")
bioRiskPlot(inputFile="risk.TCGA.txt",riskScoreFile="tcga.riskScore.pdf",survStatFile="tcga.survStat.pdf",heatmapFile="tcga.heatmap.pdf")
bioRiskPlot(inputFile="risk all.txt",riskScoreFile="all.riskScore.pdf",survStatFile="all.survStat.pdf",heatmapFile="all.heatmap.pdf")



#Fig.4 D-E

library(survival)
risk=read.table("risk.TCGA.txt",header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table("tcgaClinical.txt",sep="\t",check.names=F,header=T,row.names=1)     
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])


uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	 coxSummary = summary(cox)
	 uniTab=rbind(uniTab,
	              cbind(id=i,
	              HR=coxSummary$conf.int[,"exp(coef)"],
	              HR.95L=coxSummary$conf.int[,"lower .95"],
	              HR.95H=coxSummary$conf.int[,"upper .95"],
	              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
	              )
}
write.table(uniTab,file="tcga.uniCox.txt",sep="\t",row.names=F,quote=F)


multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="tcga.multiCox.txt",sep="\t",row.names=F,quote=F)



bioForest=function(coxFile=null,forestFile=null,forestCol=null){
		
		rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
		gene <- rownames(rt)
		hr <- sprintf("%.3f",rt$"HR")
		hrLow  <- sprintf("%.3f",rt$"HR.95L")
		hrHigh <- sprintf("%.3f",rt$"HR.95H")
		Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
		pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
		
		pdf(file=forestFile, width = 6.3,height = 4.5)
		n <- nrow(rt)
		nRow <- n+1
		ylim <- c(1,nRow)
		layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
		
		xlim = c(0,3)
		par(mar=c(4,2.5,2,1))
		plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
		text.cex=0.8
		text(0,n:1,gene,adj=0,cex=text.cex)
		text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
		text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
		
		
		par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
		xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
		plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
		arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
		abline(v=1,col="black",lty=2,lwd=2)
		boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
		points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
		axis(1)
		dev.off()
}

bioForest(coxFile="tcga.uniCox.txt",forestFile="tcga.uniForest.pdf",forestCol="green")
bioForest(coxFile="tcga.multiCox.txt",forestFile="tcga.multiForest.pdf",forestCol="red")



#Fig.4 F
library(survival)
library(regplot)
riskFile="risk all.txt"       
cliFile="all clinical for nomo.txt"   


risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "risk")]


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)


samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)


res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1<-regplot(res.cox,
              plots = c("density", "boxes"),
              clickable=F,
              title="",
              points=TRUE,
              droplines=TRUE,
              observation=rt[1,],
              rank="sd",
              failtime = c(1,3,5),
              prfail = T)