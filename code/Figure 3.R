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
source("main.R")

##############################################################
# Figure 3
##############################################################
nrepeat = 20
nsamples = 100

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

############## proportion  ##################
load("CCLE.exp.RData")
lung_exp<- RGE[,c(grep("LUNG",colnames(RGE)))]  #lung cell lines, 91 cell lines
lung_exp <- quantile_normalisation(lung_exp)    # quantile normalization

################# Fig3a
set.seed(12)
W = as.matrix(lung_exp[,1:6])
res_H = c()
res_w = c()
for (lnoise in c(0.1,0.2,0.3,0.4,0.5)){
  out = simuAll(W = W,W1.index = 1:4,csd = lnoise,nSample = nsamples,type = "GE",fsmethod = "cv",method = "all",nrep = nrepeat)
  res1 = cbind(out$RBout[,"H_AbsBias"],"Imfit",lnoise)
  res2 = cbind(out$CBSout[,"H_AbsBias"],"CBS",lnoise)
  res3 = cbind(out$RFout[,"H_AbsBias"],"RF",lnoise)
  res4 = cbind(out$PRout[,"H_AbsBias"],"PREDE",lnoise)
  res_H = rbind(res_H,res1,res2,res3,res4)
  ####
  res5 = cbind(out$RFout[,"W_corr"],"RF",lnoise)
  res6 = cbind(out$PRout[,"W_corr"],"PREDE",lnoise)
  res_w = rbind(res_w,res5,res6)
}
save(res_H,res_w,file="BiasNoise20200303.Rdata")
load("BiasNoise20200303.Rdata")
bias = as.numeric(res_H[,1])
method = as.factor(res_H[,2])
noise = as.numeric(res_H[,3])
out = data.frame(bias = bias,method = method,noise = noise)

out0 = aggregate(out$bias,by = list(out$method,factor(out$noise)),FUN = "mean")
out_sd = aggregate(out$bias,by = list(out$method,factor(out$noise)),FUN = "sd")
dat = cbind(out0,out_sd$x)
colnames(dat)= c("method","noise","bias","sd")
dat$method = factor(dat$method, levels=c("Imfit","CBS","RF", "PREDE"))
cat("calculating: p2")
p1 = ggplot(dat, aes(x=noise, y=bias, group=method,color = method)) +
  geom_line(aes(linetype=method))+
  scale_color_manual(values=brewer.pal(8, "Dark2")[c(4,3,2,1)]) +
  theme_classic()+
  geom_point(aes(shape=method)) +
  geom_errorbar(aes(ymin=bias-sd, ymax=bias+sd), width=.1) +
  xlab("Level of noise")+
  ylab("Mean absolute error")+
  labs(fill = "Method")
  #theme(legend.position = "none")

################# Fig3b
load("BiasNoise20200303.Rdata")
corr = as.numeric(res_w[,1])
method = as.factor(res_w[,2])
noise = as.numeric(res_w[,3])
out = data.frame(corr = corr,method = method,noise = noise)

out0 = aggregate(out$corr,by = list(out$method,factor(out$noise)),FUN = "mean")
out_sd = aggregate(out$corr,by = list(out$method,factor(out$noise)),FUN = "sd")

dat = cbind(out0,out_sd$x)
colnames(dat)= c("method","noise","corr","sd")
dat$method = factor(dat$method, levels=c("RF", "PREDE"))
cat("calculating: p5")
p3 = ggplot(dat, aes(x=noise, y=corr, fill=method)) +
  #geom_line(aes(linetype=method,color=method))+
  #geom_point(aes(shape=method,color=method)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=corr, ymax=corr+sd), width=.3,position=position_dodge(0.9)) +
  scale_color_manual(values=brewer.pal(8, "Dark2")[c(2,1)]) +
  theme_classic()+
  xlab("Level of noise")+
  ylab("Correlation with truth")+
  labs(fill = "Method")

################# Fig3c
set.seed(12)
W = as.matrix(lung_exp[,1:10])
res = c()
for (k1 in seq(9,5,-1)){
  out = simuAll(W = W,W1.index = 1:k1,csd = 0.1,nSample = nsamples,type = "GE",fsmethod = "cv",method = "all",nrep = nrepeat)
  res1 = cbind(out$RBout[,"H_AbsBias"],"Imfit",(10-k1)/10)
  res2 = cbind(out$CBSout[,"H_AbsBias"],"CBS",(10-k1)/10)
  res3 = cbind(out$RFout[,"H_AbsBias"],"RF",(10-k1)/10)
  res4 = cbind(out$PRout[,"H_AbsBias"],"PREDE",(10-k1)/10)
  result =rbind(res1,res2,res3,res4)
  res = rbind(res,result)
}
#save(res,file="Biasfract20200303.Rdata")
load("Biasfract20200303.Rdata")
bias = as.numeric(res[,1])
method = as.factor(res[,2])
unfract = as.numeric(res[,3])
out = data.frame(bias = bias,method = method,unfract = unfract)

out0 = aggregate(out$bias,by = list(out$method,factor(out$unfract)),FUN = "mean")
out_sd = aggregate(out$bias,by = list(out$method,out$unfract),FUN = "sd")
dat = cbind(out0,out_sd$x)
colnames(dat)= c("method","unfract","bias","sd")

dat$method = factor(dat$method, levels=c("Imfit","CBS","RF", "PREDE"))
cat("\ncalculating: p2\n")
p2 = ggplot(dat, aes(x=unfract, y=bias, group=method,color=method)) +
  geom_line(aes(linetype=method))+
  ylim(0,NA)+
  scale_color_manual(values=brewer.pal(8, "Dark2")[c(4,3,2,1)]) +
  theme_classic()+
  geom_point(aes(shape=method,color=method)) +
  geom_errorbar(aes(ymin=bias-sd, ymax=bias+sd), width=.1) +
  xlab("Proportion of unknown cell types")+
  ylab("Mean absolute error") + 
  labs(fill = "Method")

################ Fig3d
set.seed(12)
W = as.matrix(lung_exp[,1:10])

res1 = c()
for (k1 in seq(9,1,-2)){
  out = simuAll(W = W,W1.index = 1:k1,csd = 0.1,nSample = nsamples,fsmethod = "cv",method = "PR",nrep = nrepeat)
  result = cbind(out$PRout[,"W_corr"],"PREDE",(10-k1)/10)
  res1 = rbind(res1,result)
}
res2 = c()
for (k1 in seq(9,1,-2)){
  out = simuAll(W = W,W1.index = 1:k1,csd = 0.1,nSample = nsamples,fsmethod = "cv",method = "RF",nrep = nrepeat)
  result2 = cbind(out$RFout[,"W_corr"],"RF",(10-k1)/10)
  res2 = rbind(res2,result2)
}
res =rbind(res1,res2)
#save(res,file="Corfract20200303.Rdata")
load("Corfract20200303.Rdata")
cor = as.numeric(res[,1])
method = as.factor(res[,2])
unfract = as.factor(res[,3])
out = data.frame(cor = cor,method = method,unfract = unfract)

out0 = aggregate(out$cor,by = list(out$method,out$unfract),FUN = "mean")
out_sd = aggregate(out$cor,by = list(out$method,out$unfract),FUN = "sd")
                      
dat = cbind(out0,out_sd$x)
colnames(dat)= c("method","unfract","cor","sd")

dat$method = factor(dat$method, levels=c("RF", "PREDE"))
cat("\ncalculating: p4\n")
p4 = ggplot(dat, aes(x=unfract, y=cor, fill=method)) +
  #geom_line(aes(linetype=method,color=method))+
  #geom_point(aes(shape=method,color=method)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  #ylim(0.5,1)+
  theme_classic()+
  geom_errorbar(aes(ymin=cor, ymax=cor+sd), width=.3,position=position_dodge(0.9)) +
  scale_color_manual(values=brewer.pal(8, "Dark2")[c(2,1)]) +
  xlab("Proportion of unknown cell types")+
  ylab("Correlation with truth") +
  labs(fill = "Method")


p =ggarrange(p1,p2,p3,p4,ncol=2,nrow=2,labels=c("A","B","C","D"))
ggsave(plot = p,filename="Figure2.pdf",width=7, height=4.3)
