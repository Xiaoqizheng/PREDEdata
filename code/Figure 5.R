rm(list = ls())
setwd("/Users/wwzhang/Desktop/PREDE_Figure5/")
load("Figure5.RData")
#*************************
## missing one cell type
#************************
## correlation matrix of predicted expression vs. true expression of missing cell type
expr_data1 = matrix(NA,nrow = 3,ncol=3)
colnames(expr_data1) = c("Liver","Brain","Lung") ## missing cell type
rownames(expr_data1) = c("PREDE","DeMixT","ISOpure") ## methods

## ISOpure
expr_data1[3,] = diag(cor(ISOpure.out$pf,pure_mean,method = "pearson")) 
prop_ISOpure = ISOpure.out$pp

## PREDE
exprs = cbind(out.liver.PR$W[,1]*sc,out.brain.PR$W[,2]*sc,out.lung.PR$W[,3]*sc)
prop_PREDE_liver = t(out.liver.PR$H)
prop_PREDE_brain = t(out.brain.PR$H)
prop_PREDE_lung = t(out.lung.PR$H)
expr_data1[1,] = diag(cor(exprs,pure_mean,method = "pearson"))

## DeMixT
expr_data1[2,] = diag(cor(cbind(apply(out.liver$ExprT,1,mean),apply(out.brain$ExprT,1,mean),
                                apply(out.lung$ExprT,1,mean)),pure_mean,method="pearson")) 
prop_DeMixT_liver = t(rbind(1-apply(out.liver$pi,2,sum),out.liver$pi))
prop_DeMixT_brain = t(rbind(out.brain$pi[1,],1-apply(out.brain$pi,2,sum),out.brain$pi[2,]))
prop_DeMixT_lung = t(rbind(out.lung$pi,1-apply(out.lung$pi,2,sum)))


## plot
library(RColorBrewer)
color = c(brewer.pal(9, "Set1")[1],"dodgerblue4",brewer.pal(9, "Set1")[3])

pdf("Figure5.pdf",width = 10,height = 5)
par(mar = c(4, 3, 2, 1.1), mgp = c(1.9, 0.5, 0),mfrow = c(2,4))

plot(trueProp[,1],prop_PREDE_liver[,1],xlim = c(0,1),ylim = c(0,1),main = "Liver unknown",
     cex = 0.5,col = color[1],pch = 19,xlab = "True Proportion", ylab = "Estimated Proportion",
     cex.main = 1,cex.lab=1,cex.axis=1,frame.plot=F)
points(trueProp[,2],prop_PREDE_liver[,2],col = color[1],cex = 0.5,pch = 3)
points(trueProp[,3],prop_PREDE_liver[,3],col = color[1],cex = 0.5,pch = 17)

points(trueProp[,1],prop_DeMixT_liver[,1],col = color[2],cex = 0.5,pch = 19)
points(trueProp[,2],prop_DeMixT_liver[,2],col = color[2],cex = 0.5,pch = 3)
points(trueProp[,3],prop_DeMixT_liver[,3],col = color[2],cex = 0.5,pch = 17)

points(trueProp[,1],prop_ISOpure[,1],col = color[3],cex = 0.5,pch = 19)

abline(0,1,col = "gray40",lwd = 1.5,lty=1) 
legend("bottomright",legend=c("Liver","Brain","Lung"),pch = c(19,3,17), bty="n",cex=1)

MAE_PREDE = round(mean(colSums(abs(prop_PREDE_liver - trueProp)) / nrow(trueProp)),3)
MAE_DeMixT = round(mean(colSums(abs(prop_DeMixT_liver - trueProp)) / nrow(trueProp)),3)
MAE_ISOpure = round(sum(abs(prop_ISOpure[,1]-trueProp[,1]))/nrow(trueProp),3)
text(paste0("PREDE(MAE=",MAE_PREDE,")"),x = 0.28,y = 1,cex = 1,col=color[1])
text(paste0("DeMixT(MAE=",MAE_DeMixT,")"),x = 0.28,y = 0.9,cex = 1,col=color[2])
text(paste0("ISOpure(MAE=",MAE_ISOpure,")"),x = 0.28,y = 0.8,cex = 1,col=color[3])

## brain
plot(trueProp[,1],prop_PREDE_brain[,1],xlim = c(0,1),ylim = c(0,1),main = "Brain unknown",
     cex = 0.5,col = color[1],pch = 19,xlab = "True Proportion", ylab = "Estimated Proportion",
     cex.main = 1,cex.lab=1,cex.axis=1,frame.plot=F)
points(trueProp[,2],prop_PREDE_brain[,2],col = color[1],cex = 0.5,pch = 3)
points(trueProp[,3],prop_PREDE_brain[,3],col = color[1],cex = 0.5,pch = 17)

points(trueProp[,1],prop_DeMixT_brain[,1],col = color[2],cex = 0.5,pch = 19)
points(trueProp[,2],prop_DeMixT_brain[,2],col = color[2],cex = 0.5,pch = 3)
points(trueProp[,3],prop_DeMixT_brain[,3],col = color[2],cex = 0.5,pch = 17)

points(trueProp[,2],prop_ISOpure[,2],col = color[3],cex = 0.5,pch = 3)

abline(0,1,col = "gray40",lwd = 1.5,lty=1) 
legend("bottomright",legend=c("Liver","Brain","Lung"),pch = c(19,3,17), bty="n",cex=1)

MAE_PREDE = round(mean(colSums(abs(prop_PREDE_brain - trueProp)) / nrow(trueProp)),3)
MAE_DeMixT = round(mean(colSums(abs(prop_DeMixT_brain - trueProp)) / nrow(trueProp)),3)
MAE_ISOpure = round(sum(abs(prop_ISOpure[,2]-trueProp[,2]))/nrow(trueProp),3)
text(paste0("PREDE(MAE=",MAE_PREDE,")"),x = 0.28,y = 1,cex = 1,col=color[1])
text(paste0("DeMixT(MAE=",MAE_DeMixT,")"),x = 0.28,y = 0.9,cex = 1,col=color[2])
text(paste0("ISOpure(MAE=",MAE_ISOpure,")"),x = 0.28,y = 0.8,cex = 1,col=color[3])

## lung
plot(trueProp[,1],prop_PREDE_lung[,1],xlim = c(0,1),ylim = c(0,1),main = "Lung unknown",
     cex = 0.5,col = color[1],pch = 19,xlab = "True Proportion", ylab = "Estimated Proportion",
     cex.main = 1,cex.lab=1,cex.axis=1,frame.plot=F)
points(trueProp[,2],prop_PREDE_lung[,2],col = color[1],cex = 0.5,pch = 3)
points(trueProp[,3],prop_PREDE_lung[,3],col = color[1],cex = 0.5,pch = 17)

points(trueProp[,1],prop_DeMixT_lung[,1],col = color[2],cex = 0.5,pch = 19)
points(trueProp[,2],prop_DeMixT_lung[,2],col = color[2],cex = 0.5,pch = 3)
points(trueProp[,3],prop_DeMixT_lung[,3],col = color[2],cex = 0.5,pch = 17)

points(trueProp[,3],prop_ISOpure[,3],col = color[3],cex = 0.5,pch = 17)

abline(0,1,col = "gray40",lwd = 1.5,lty=1) 
legend("bottomright",legend=c("Liver","Brain","Lung"),pch = c(19,3,17), bty="n",cex=1)


MAE_PREDE = round(mean(colSums(abs(prop_PREDE_lung - trueProp)) / nrow(trueProp)),3)
MAE_DeMixT = round(mean(colSums(abs(prop_DeMixT_lung - trueProp)) / nrow(trueProp)),3)
MAE_ISOpure = round(sum(abs(prop_ISOpure[,3]-trueProp[,3]))/nrow(trueProp),3)
text(paste0("PREDE(MAE=",MAE_PREDE,")"),x = 0.28,y = 1,cex = 1,col=color[1])
text(paste0("DeMixT(MAE=",MAE_DeMixT,")"),x = 0.28,y = 0.9,cex = 1,col=color[2])
text(paste0("ISOpure(MAE=",MAE_ISOpure,")"),x = 0.28,y = 0.8,cex = 1,col=color[3])

## barplot of correltion between expression
x <- barplot(expr_data1,beside=TRUE,space=c(0,2),col=color,cex.axis=1,xaxt="n",las=1,
             ylab="Correlation of expression",cex.lab=1,ylim=c(0,1),xlab="")
axis(1,x[2,],cex.axis=1,labels = F ,tck=-0.02)
text(cex=1, x=x[2,], y=-0.1, colnames(expr_data1), xpd=TRUE,srt=0)
legend(6,1.17,rownames(expr_data1),fill=color,bty="n",cex=1)

#*********************************
## missing two cell types
#*********************************
expr_data2 = matrix(NA,nrow=2,ncol=3)
expr_data2[,1] = diag(cor(out.brain.lung.PR$W*sc,pure_mean))[c(2,3)]
prop_PREDE_brain_lung = t(out.brain.lung.PR$H)
expr_data2[,2] = diag(cor(out.liver.brain.PR$W*sc,pure_mean))[c(1,2)]
prop_PREDE_liver_brain = t(out.liver.brain.PR$H)
expr_data2[,3] = diag(cor(out.liver.lung.PR$W*sc,pure_mean))[c(1,3)]
prop_PREDE_liver_lung = t(out.liver.lung.PR$H)
colnames(expr_data2) = c("Brain_Lung","Liver_Brain","Liver_Lung")

## plot
## Brain and Lung known
plot(trueProp[,1],prop_PREDE_brain_lung[,1],xlim = c(0,1),ylim = c(0,1),main = "Brain & Lung unknown",
     cex = 0.5,col = color[1],pch = 19,xlab = "True Proportion", ylab = "Estimated Proportion",
     cex.main = 1,cex.lab=1,cex.axis=1,frame.plot = F)
points(trueProp[,2],prop_PREDE_brain_lung[,2],col = color[1],cex = 0.5,pch = 3)
points(trueProp[,3],prop_PREDE_brain_lung[,3],col = color[1],cex = 0.5,pch = 17)
abline(0,1,col = "gray40",lwd = 1.5,lty=1) 
legend("bottomright",legend=c("Liver","Brain","Lung"),pch = c(19,3,17), bty="n",cex=1)
MAE_PREDE = round(mean(colSums(abs(prop_PREDE_brain_lung - trueProp)) / nrow(trueProp)),3)
text(paste0("PREDE(MAE=",MAE_PREDE,")"),x = 0.28,y = 0.9,cex = 1,col=color[1])

## Liver and Brain known
plot(trueProp[,1],prop_PREDE_liver_brain[,1],xlim = c(0,1),ylim = c(0,1),main = "Liver & Brain unknown",
     cex = 0.5,col = color[1],pch = 19,xlab = "True Proportion", ylab = "Estimated Proportion",
     cex.main = 1,cex.lab=1,cex.axis=1,frame.plot=F)
points(trueProp[,2],prop_PREDE_liver_brain[,2],col = color[1],cex = 0.5,pch = 3)
points(trueProp[,3],prop_PREDE_liver_brain[,3],col = color[1],cex = 0.5,pch = 17)
abline(0,1,col = "gray40",lwd = 1.5,lty=1) 
legend("bottomright",legend=c("Liver","Brain","Lung"),pch = c(19,3,17), bty="n",cex=1)
MAE_PREDE = round(mean(colSums(abs(prop_PREDE_liver_brain - trueProp)) / nrow(trueProp)),3)
text(paste0("PREDE(MAE=",MAE_PREDE,")"),x = 0.28,y = 0.9,cex = 1,col=color[1])

## Liver and Lung known
plot(trueProp[,1],prop_PREDE_liver_lung[,1],xlim = c(0,1),ylim = c(0,1),main = "Liver & lung unknown",
     cex = 0.5,col = color[1],pch = 19,xlab = "True Proportion", ylab = "Estimated Proportion",
     cex.main = 1,cex.lab=1,cex.axis=1,frame.plot = F)
points(trueProp[,2],prop_PREDE_liver_lung[,2],col = color[1],cex = 0.5,pch = 3)
points(trueProp[,3],prop_PREDE_liver_lung[,3],col = color[1],cex = 0.5,pch = 17)
abline(0,1,col = "gray40",lwd = 1.5,lty=1) 
legend("bottomright",legend=c("Liver","Brain","Lung"),pch = c(19,3,17), bty="n",cex=1)
MAE_PREDE = round(mean(colSums(abs(prop_PREDE_liver_lung - trueProp)) / nrow(trueProp)),3)
text(paste0("PREDE(MAE=",MAE_PREDE,")"),x = 0.28,y = 0.9,cex = 1,col=color[1])

color2 = c(brewer.pal(9, "RdGy")[c(1,3,4)])
x <- barplot(expr_data2,beside=TRUE,space=c(0,2),col=c(color2[2],color2[3],color2[1],color2[2],color2[1],color2[3]),
             cex.axis=1,xaxt="n",las=1,ylab="Correlation of expression",cex.lab=1,ylim=c(0,1),xlab="")
axis(1,x[2,]-0.5,cex.axis=1,labels = F ,tck=-0.02)
text(cex=1, x=x[2,]-0.5, y=-0.1, colnames(expr_data2), xpd=TRUE,srt=60)
legend(1.3,1.28,c("Liver","Brain","Lung"),fill = color2, bty="n",cex=1)
dev.off()









