library(survival)
library(survminer)
library(forestplot)


rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)       
rt$futime=rt$futime/365
gene=colnames(rt)[3]
pFilter=0.05            


for(i in levels(factor(rt[,"CancerType"]))){
	rt1=rt[(rt[,"CancerType"]==i),]
	res.cut=surv_cutpoint(rt1, time = "futime", event = "fustat", variables =gene)
	res.cat=surv_categorize(res.cut)
	#fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
	#group=ifelse(rt1[,gene]>res.cut,"high","low")
	diff=survdiff(Surv(futime, fustat) ~riskScore,data = res.cat)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<pFilter){
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		fit <- survfit(Surv(futime, fustat) ~ riskScore, data = res.cat)
		
		surPlot=ggsurvplot(fit, 
				    data=res.cat,
				    title=paste0("Cancer: ",i),
				    pval=pValue,
				    pval.size=6,
				    legend.labs=c("high","low"),
				    legend.title=paste0(gene," levels"),
				    font.legend=12,
				    xlab="Time(years)",
				    ylab="Overall survival",
				    break.time.by = 1,
				    palette=c("red","blue"),
				    conf.int=F,
				    fontsize=4,
				    risk.table=TRUE,
				    risk.table.title="",
				    risk.table.height=.25)
		pdf(file=paste0("OS survival.",i,".pdf"),onefile = FALSE,
				    width = 6,             
				    height =5)             #
		print(surPlot)
		dev.off()
	}
}


rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)    
rt$futime=rt$futime/365
gene=colnames(rt)[3]


outTab=data.frame()
for(i in levels(factor(rt[,"CancerType"]))){
	rt1=rt[(rt[,"CancerType"]==i),]
	cox=coxph(Surv(futime, fustat) ~ rt1[,gene], data = rt1)
	coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	outTab=rbind(outTab,
	             cbind(cancer=i,
	                   HR=coxSummary$conf.int[,"exp(coef)"],
	                   HR.95L=coxSummary$conf.int[,"lower .95"],
	                   HR.95H=coxSummary$conf.int[,"upper .95"],
			           pvalue=coxP) )
}
write.table(outTab,file="OS cox.result.txt",sep="\t",row.names=F,quote=F)    



bioForest=function(coxFile=null, forestFile=null){
	rt <- read.table(coxFile,header=T,sep="\t",check.names=F,row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrLow[hrLow<0.001]=0.001
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	
	height=(nrow(rt)/15)+5
	pdf(file=forestFile, width=8, height=height)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
		
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	LOGindex=10 
	hrLow = log(as.numeric(hrLow),LOGindex)
	hrHigh = log(as.numeric(hrHigh),LOGindex)
	hr = log(as.numeric(hr),LOGindex)
	xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), "red", "blue")
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
	a1 = axis(1,labels=F,tick=F)
	axis(1,a1,LOGindex^a1)
	dev.off()
}

bioForest(coxFile="OS cox.result.txt", forestFile="OS forest.pdf")