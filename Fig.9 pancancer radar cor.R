library(fmsb) 

data=read.table("TMB fmsbInput.txt",header=T,sep="\t",row.names=1,check.names=F)   
maxValue=ceiling(max(abs(data))*10)/10
data=rbind(rep(maxValue,ncol(data)),rep(-maxValue,ncol(data)),data)

colors="red"
corStat=read.table("TMB corStat.txt",header=T,sep="\t",row.names=1,check.names=F)
colnames(data)=paste0(colnames(data),corStat$sig)


pdf(file="TMB radar.pdf",height=7,width=7)
radarchart( data, axistype=1 , 
    pcol=colors,                 
    plwd=2 ,                     
    plty=1,                      
    cglcol="grey",               
    cglty=1,                      
    caxislabels=seq(-maxValue,maxValue,maxValue/2),    
    cglwd=1.2,                 
    axislabcol="blue",           
    vlcex=0.8                    
)
dev.off()