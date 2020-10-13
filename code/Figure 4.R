rm(list = ls())
library(quadprog)
library(e1071)
library(corrplot)
library(RColorBrewer)
library(gplots)
library(ggpubr)
library(vioplot)
library(ggplot2)
library(preprocessCore)
source("functions.R")

##############################################################
# Figure 4
##############################################################
pdf("Figure4.pdf",width = 9.7,height = 5)
par(mar = c(3.5, 3, 1.6, 1.1), mgp = c(1.9, 0.5, 0),mfrow = c(2,4))

########## read gene expression profile of mixed cells
load("Mix_cell.Rdata")
cellName = colnames( Mix_cell)

##########adjust the proportion of cell lines
B_ca = generate_bulk(celllines = W[,1:3],nSample = 100,csd = 0.1)
B_im = generate_bulk(celllines = W[,4:6],nSample = 100,csd = 0.1)
B_nm = generate_bulk(celllines = W[,7:8],nSample = 100,csd = 0.1)

B_Y = B_ca$Y*0.6+B_im$Y*0.2+B_nm$Y*0.2             # cancer 60%,immune 20%,nornal 20%
B_H = rbind(B_ca$H*0.6,B_im$H*0.2,B_nm$H*0.2)

########## choose "cv" to select features

rownames(B_Y) <- rownames(Mix_cell)
rownames(W) <- rownames(Mix_cell)
difgene = select_feature(B_Y,method = "cv",nmarker=1000, startn=0)
                                                   # use "cv" to select features
################unknown cancer cell
out.PR <- PR(Y = B_Y[difgene,], W = W[difgene,], W1 = W[difgene,c(4:6,7:8)],K = ncol(W),type = "GE",iters = 500,rssDiffStop=1e-6)
for(i in 1:3){
    plot(out.PR$W[,i],W[difgene,i],xlim = c(0,15),ylim = c(0,15),col = "#00000050",pch = 19,cex.main=1,xlab = "Estimated expression profile",ylab = "True expression profile",main = cellName[i])
    text(2,14.5,labels = paste0("R = ",round(cor(out.PR$W[,i],W[difgene,i]),3)))
    abline(0,1,lty = 2,lwd = 2) 
}

plot(out.PR$H[1,],B_H[1,],xlim = c(0,0.5),ylim = c(0,0.5),cex = 0.6,col = "blue",pch = 19,xlab = "Estimated proportion",ylab = "True proportion",)
points(out.PR$H[2,],B_H[2,],col = adjustcolor("red", alpha.f = 0.8),cex = 0.6,pch = 3)
points(out.PR$H[3,],B_H[3,],col = adjustcolor( "green", alpha.f = 0.8),cex = 0.6,pch = 17)
abline(0,1,lty = 2,col = "gray40",lwd = 1.5) 
legend("topleft",legend=c("AU565","BT474","BT549"), col=c("blue","red","green"),
       pch=c(19,3,17), ncol=1,bty="n",cex=0.9)

################unknown immune cell
out.PR.im <- PR(Y = B_Y[difgene,], W = W[difgene,], W1 = W[difgene,c(1:3,7:8)],K = ncol(W),type = "GE",iters = 500,rssDiffStop=1e-6)
for(i in 4:6){
    plot(out.PR.im$W[,i],W[difgene,i],xlim = c(0,15),ylim = c(0,15),col = "#00000050",cex.main=1,pch = 19,xlab = "Estimated expression profile",ylab = "True expression profile",main.cex =0.9,main = cellName[i])
    text(2,14.5,labels = paste0("R = ",round(cor(out.PR.im$W[,i],W[difgene,i]),3)))
    abline(0,1,lty = 2,lwd = 2) 
}

plot(out.PR$H[4,],B_H[4,],xlim = c(0,0.2),ylim = c(0,0.2),cex = 0.6,col = "blue",pch = 19,xlab = "Estimated proportion",ylab = "True proportion")
points(out.PR$H[5,],B_H[5,],col = adjustcolor( "red", alpha.f = 0.8),cex = 0.6,pch = 3)
points(out.PR$H[6,],B_H[6,],col = adjustcolor( "green", alpha.f = 0.8),pch = 17,cex = 0.6)
abline(0,1,lty = 2,col = "gray40",lwd = 1.5) 
legend("topleft",legend=c("CD8 T cell","CD4 T cell","B cell"), col=c("blue","red","green"),
       pch=c(19,3,17), ncol=1,bty="n",cex=0.9)

dev.off()




