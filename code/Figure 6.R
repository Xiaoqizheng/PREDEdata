library(RColorBrewer)
load("Figure6.RData")

################################################################
## plot barplot
colors = brewer.pal(8, "Dark2")[c(4,7,3,2,1)]
pdf("Figure6.pdf", height = 2.5*2, width = 2.5*3)
par(mar = c(3.6, 3.5, 2.5, 1.1), mgp = c(2, 0.5, 0), mfrow = c(2,3))

x <- barplot(data_known_result[[1]][c(1,5,2,3,4),], beside=TRUE,col = colors,cex.axis = 1,xaxt="n",las=1, xlab="", 
             ylab="Pearson Correlation",cex.lab = 1,space=c(0,2),border=NA,ylim = c(min(data_known_result[[1]]),1),
             main = "3 cell types known",cex.main = 1)
text(cex = 1, x=x[3,]-0.5, y=-0.14,colnames(data_known_result[[1]][c(1,5,2,3,4),]), xpd=TRUE,srt = 90)
legend(8,1.15,rownames(data_known_result[[1]][c(1,5,2,3,4),]),fill = colors,bty = "n",cex = 1)

x <- barplot(data_known_result[[2]][c(1,5,2,3,4),], beside=TRUE,col = colors,cex.axis = 1,xaxt="n",las=1, xlab="", 
             ylab="Pearson Correlation",cex.lab = 1,space=c(0,2),border=NA,ylim = c(min(data_known_result[[2]]),1),
             main = "5 cell types known",cex.main = 1)
text(cex = 1, x=x[2,]-0.5, y=-0.14,colnames(data_known_result[[2]][c(1,5,2,3,4),]), xpd=TRUE,srt = 90)

x <- barplot(data_known_result[[3]][c(1,5,2,3,4),], beside=TRUE,col = colors,cex.axis = 1,xaxt="n",las=1, xlab="", 
             ylab="Pearson Correlation",cex.lab = 1,space=c(0,2),border=NA,ylim = c(min(data_known_result[[3]]),1),
             main = "7 cell types known",cex.main = 1)
text(cex = 1, x=x[2,]-0.5, y=-0.14,colnames(data_known_result[[3]][c(1,5,2,3,4),]), xpd=TRUE,srt = 90)

## scatter
colors = brewer.pal(8, "Dark2")[c(4,7,3,2,1,5)]
plot(trueProp1[,4],data_unknown_result[[1]][,1],cex = 0.6,col = colors[1],pch = 19,xlab = "True Proportion", 
     ylab = "PREDE Proportion",cex.main = 1,xlim = c(0,45),ylim = c(0,max(data_unknown_result[[1]])),
     frame.plot=F,main = "3 cell types known")
points(trueProp1[,5],data_unknown_result[[1]][,2],col = colors[2],cex = 0.6,pch = 19)
points(trueProp1[,6],data_unknown_result[[1]][,3],col = colors[3],cex = 0.6,pch = 19)
points(trueProp1[,7],data_unknown_result[[1]][,4],col = colors[4],cex = 0.6,pch = 19)
points(trueProp1[,8],data_unknown_result[[1]][,5],col = colors[5],cex = 0.6,pch = 19)
points(trueProp1[,9],data_unknown_result[[1]][,6],col = colors[6],cex = 0.6,pch = 19)
abline(0,1,lty = 2,col = "gray40",lwd = 1.5) 
legend("bottomright",legend=colnames(trueProp1)[4:9], col = colors,
       pch = rep(19,6), ncol=1,bty = "n",cex = 0.8)
MC = round(mean(diag(cor(trueProp1[,4:9],data_unknown_result[[1]]))),2)
text(paste0("R = ",MC),x = 15,y = 55,cex = 1)

plot(trueProp1[,6],data_unknown_result[[2]][,1],cex = 0.6,col = colors[1],pch = 19,xlab = "True Proportion", 
     ylab = "PREDE Proportion",cex.main = 1,xlim = c(0,45),ylim = c(0,max(data_unknown_result[[1]])),
     frame.plot=F,main = "5 cell types known")
points(trueProp1[,7],data_unknown_result[[2]][,2],col = colors[2],cex = 0.6,pch = 19)
points(trueProp1[,8],data_unknown_result[[2]][,3],col = colors[3],cex = 0.6,pch = 19)
points(trueProp1[,9],data_unknown_result[[2]][,4],col = colors[4],cex = 0.6,pch = 19)
abline(0,1,lty = 2,col = "gray40",lwd = 1.5) 
legend("bottomright",legend=colnames(trueProp1)[6:9], col = colors[1:4],
       pch = rep(19,4), ncol=1,bty = "n",cex = 0.8)
MC = round(mean(diag(cor(trueProp1[,6:9],data_unknown_result[[2]]))),2)
text(paste0("R = ",MC),x = 15,y = 55,cex = 1)


plot(trueProp1[,8],data_unknown_result[[3]][,1],cex = 0.6,col = colors[1],pch = 19,xlab = "True Proportion", 
     ylab = "PREDE Proportion",cex.main = 1,xlim = c(0,45),ylim = c(0,max(data_unknown_result[[1]])),
     frame.plot=F,main = "7 cell types known")
points(trueProp1[,9],data_unknown_result[[3]][,2],col = colors[2],cex = 0.6,pch = 19)
abline(0,1,lty = 2,col = "gray40",lwd = 1.5) 
legend("bottomright",legend=colnames(trueProp1)[8:9], col = colors[1:2],
       pch = rep(19,2), ncol=1,bty = "n",cex = 0.8)
MC = round(mean(diag(cor(trueProp1[,8:9],data_unknown_result[[3]]))),2)
text(paste0("R = ",MC),x = 15,y = 55,cex = 1)
dev.off()







