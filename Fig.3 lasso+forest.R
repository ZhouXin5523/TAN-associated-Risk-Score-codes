#Fig.3A


library("glmnet")
library("survival")

#LASSO
rt=read.table("uniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)       

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit<-glmnet(x,y,family = "cox", maxit = 1000)
plot(fit,xvar="lambda",label=T)
dev.off()

fitcv <- cv.glmnet(x,y,family="cox", alpha=1,nfolds=10)
plot(fitcv)
coef(fitcv, s="lambda.min")

coef.min = coef(fitcv, s = "lambda.min")
active.min = which(coef.min != 0)
index.min = coef.min[active.min]
row.names(coef.min)[active.min]

model <- coxph(s, data = rt)
summary(model,data=rt)

RiskScore<-predict(model,type = "risk")
names(RiskScore) = rownames(rt)

dev.off()
#Fig.3B

ggforest(model,  
         data =rt,  
         main = 'Hazard ratio ', 
         cpositions = c(0.05, 0.15, 0.35),  
         fontsize = 1, 
         refLabel = 'reference', 
         noDigits = 3 
)

dev.off()
