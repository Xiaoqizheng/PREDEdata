rm(list = ls())
library(RColorBrewer)
library(gplots)
library(factoextra)     # A package which computes 26 distance
library(vegan)          # A package which computes Bray-curtis distance
library(quadprog)
library(splines)
library(survival)
library(GGally)
##########################################
#####Figure7A barplot 
##########################################
load("Figure7a.RData")
load("Figure7b.RData")
load("Figure7c.RData")

rownames(data.BRCA.barplot) <-c("B cell","CD4 T cell","CD8 T cell",
                                 "NK","Neutrophil","Macrophage","Dendritic cell")
colnames(data.BLCA.barplot) <- c("Infiltrated","Squamous","Papillary","Neuronal","Luminal")
colnames(data.SKCM.barplot) <- c("Type3","Type2","Type1","Type4")
data.bar = cbind(data.BRCA.barplot,data.SKCM.barplot,data.BLCA.barplot)
for (i in 1:ncol(data.bar)){
    data.bar[,i] = (data.bar[,i]/colSums(data.bar)[i])*100
}

colors = c(brewer.pal(9, "Blues")[8],brewer.pal(9, "Greens")[5],brewer.pal(9, "Set3")[5],
           brewer.pal(9, "YlOrRd")[4],brewer.pal(9, "Set1")[1],brewer.pal(9, "Purples")[6],
           brewer.pal(8, "Set2")[7])

pdf("Figure 7a.pdf",width=5,heigh=4)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
x<-barplot(data.bar,col = colors,ylim =  c(0,100),space = c(rep(0,5),0.8,rep(0,3),0.8,rep(0,4)),
        xlab = "",ylab="Estimated proportion of immune cells(%)",las = 2,cex.lab = 1,
        cex.names = 0.9,cex.axis = 0.8)
legend("right", inset=-0.48, legend=rownames(data.bar),bty="n",fill=colors)
axis(1,x,cex.axis = 0.8,labels = F ,tck=-0.02)
text(cex = 1, x=c(2.6, 8, 13), y =100,labels = c("BRCA", "SKCM","BLCA"),pos=3,xpd = NA)
dev.off()

##########################################
####Figure7B heatmap
##########################################
load("Figure7a.RData")

pdf("Figure7b.pdf",width=5.5,heigh=5)
colnames(data.BRCA.heatmap.samp) <- c("B cell","CD4 T cell","CD8 T cell","NK","Neutrophil","Macrophage","Dendritic cell",
                                      "New detected cell 1", "New detected cell 2",  "New detected cell 3",  "New detected cell 4", 
                                      "New detected cell 5", "New detected cell 6")
data.BRCA.heatmap.samp = data.BRCA.heatmap.samp[,c(8:13,1:7)]
colors    = c(brewer.pal(9, "Blues")[8],brewer.pal(9, "Reds")[6], brewer.pal(9, "Oranges")[5],
              brewer.pal(9, "Greens")[7],brewer.pal(9, "Purples")[7])
nsubty = c(508,169,190,35,78)
ColColors = c(rep(brewer.pal(9, "Blues")[8], nsubty[1]), rep(brewer.pal(9, "Reds")[6], nsubty[2]),
              rep(brewer.pal(9, "Oranges")[5], nsubty[3]),rep(brewer.pal(9, "Greens")[7], nsubty[4]),
              rep(brewer.pal(9, "Purples")[7], nsubty[5]))

heatmap.2(t(data.BRCA.heatmap.samp), scale = "none", col=brewer.pal(9, "Reds"),dendrogram="none",density.info="none",
          trace="none", ColSideColors =ColColors,labCol = FALSE,Rowv=FALSE, cexRow = 0.8, key=T, keysize = 1.5,key.xlab = "Proportion",
          distfun= function(x) vegdist(x,'bray'))
legend("top",inset=0.05,  legend=c("LumA", "Basal","LumB", "Her2", "Normal-like"),fill=colors,bty="n",ncol=5,cex=0.6)
dev.off()


##########################################
####Survival curve
##########################################

pdf("Figure 5c.pdf",width=3*3,height = 3*2)
par(mar = c(3.5, 4.5, 1.6, 1.1), mgp = c(1.9, 0.5, 0),mfcol=c(2,3))
##########################################
## load BRCA clinical data by TCGAbiolinks
load("BRCA_Subty.Rdata")

## match samples
proportions = data.BRCA.heatmap.samp
rownames(proportions) = gsub(pattern = ".",replacement = "-",x = rownames(proportions),fixed = T)
isc = intersect(substr(rownames(proportions),1,12),rownames(BRCA_Subty))
id2 = match(isc,rownames(BRCA_Subty))
id1 = match(isc,substr(rownames(proportions),1,12))
BRCA_Subty = BRCA_Subty[id2,]      #980 samples
proportions=proportions[id1,]      #980 samples

## Data for survival analysis
BRCA_Subty$days_to_death[which(BRCA_Subty$days_to_death=="NA")]<-0
BRCA_Subty$days_to_last_followup[which(BRCA_Subty$days_to_last_followup=="NA")]<-0
for (i in 1:nrow(BRCA_Subty)){
  BRCA_Subty$days_to_death[i]<-max(BRCA_Subty$days_to_death[i],BRCA_Subty$days_to_last_followup[i])
}
index = which(BRCA_Subty$days_to_death!="0") #969 samples
BRCA_Subty = BRCA_Subty[index,]          
proportions =proportions[index,]
BRCA_Subty$days_to_death = as.numeric(BRCA_Subty$days_to_death)
allData = BRCA_Subty
allData$vital_status= allData$vital_status== "Dead"

for(j in c(6,10) ){
  thisName = colnames(proportions)[j]
  allData$group = "Low"
  #allData$group[proportions[,thisName] > median(proportions[,thisName])] = "High"
  quantiles = quantile(proportions[,j], prob = seq(0, 1, length = 6))
  m1 = quantiles[2]
  m2 = quantiles[5]
  allData$group[proportions[,thisName] <= m1] = "group1"
  allData$group[proportions[,thisName] >= m2] = "group2"
  allData01 = subset(allData,group == "group1")
  allData02 = subset(allData,group == "group2") 
  allData0 = rbind(allData01,allData02)
  SurvObj = with(allData0,Surv(days_to_death, vital_status))
  fit = coxph(SurvObj ~ as.factor(group), allData0)
  pval = round(summary(fit)$coefficients[5],4)
  plot(survfit(SurvObj ~ as.factor(group), data = allData0), col=c("red", "blue"),
       mark.time=TRUE, xlab="Days", ylab="Survivial rate", main = paste0("BRCA:", thisName))
  legend("topright", legend=c("Top 20%", "Bottom 20%"),
         col=c("red", "blue"), lty=c(1,1), pch = c("+","+"),bty="n",cex=1.1)
  legend("bottomleft",legend=paste0("p=",pval), bty="n",cex=1.1)
}


##########################################
## load SKCM clinical data and the predicted proportion by PREDE
load("SKCM_Subty.Rdata")

## match samples
proportions = data.SKCM.heatmap.samp
rownames(proportions) = gsub(pattern = ".",replacement = "-",x = rownames(proportions),fixed = T)
isc = intersect(substr(rownames(proportions),1,12),rownames(SKCM_Subty))
id2 = match(isc,rownames(SKCM_Subty))
id1 = match(isc,substr(rownames(proportions),1,12))
SKCM_Subty = SKCM_Subty[id2,]     #319 samples
proportions=proportions[id1,]     #319 samples

## Data for survival analysis
index = which(SKCM_Subty$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU!="-"&
                SKCM_Subty$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU!="[Not Available]")
SKCM_Subty =SKCM_Subty[index,] 
proportions=proportions[index,]
SKCM_Subty$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU = as.character(SKCM_Subty$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU)
SKCM_Subty$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU = as.numeric(SKCM_Subty$CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU)
allData = SKCM_Subty
allData$CURATED_VITAL_STATUS= allData$CURATED_VITAL_STATUS== "Dead"

for(j in c(5,9) ){
  thisName = colnames(proportions)[j]
  allData$group = "Low"
  quantiles = quantile(proportions[,j], prob = seq(0, 1, length = 6))
  m1 = quantiles[2]
  m2 = quantiles[5]
  allData$group[proportions[,thisName] <= m1] = "group1"
  allData$group[proportions[,thisName] >= m2] = "group2"
  allData01 = subset(allData,group == "group1")
  allData02 = subset(allData,group == "group2") 
  allData0 = rbind(allData01,allData02)
  SurvObj = with(allData0, Surv(CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU, CURATED_VITAL_STATUS))
  fit = coxph(SurvObj ~ as.factor(group), allData0)
  pval = round(summary(fit)$coefficients[5],4)
  plot(survfit(SurvObj ~ as.factor(group), data = allData0), col=c("red", "blue"),
         mark.time=TRUE, xlab="Days", ylab="Survivial rate", main = paste0("SKCM:",thisName))
  legend("topright", legend=c("Top 20%", "Bottom 20%"),
         col=c("red", "blue"), lty=c(1,1), pch = c("+","+"),bty="n",cex=1.1)
  legend("bottomleft",legend=paste0("p=",pval), bty="n",cex=1.1)
}

##########################################
## load clinical data and the predicted proportion by PREDE
load("BLCA_Subty.Rdata")

## match samples
proportions = data.BLCA.heatmap.samp
rownames(proportions) = gsub(pattern = ".",replacement = "-",x = rownames(proportions),fixed = T)
isc = intersect(substr(rownames(proportions),1,12),rownames(BLCA_Subty))
id2 = match(isc,rownames(BLCA_Subty))
id1 = match(isc,substr(rownames(proportions),1,12))
BLCA_Subty = BLCA_Subty[id2,]      
proportions=proportions[id1,]

## Data for survival analysis
BLCA_Subty$`Days until death`[which(BLCA_Subty$`Days until death`=="NA")]<-0
BLCA_Subty$`Days to last followup`[which(BLCA_Subty$`Days to last followup`=="NA")]<-0
for (i in 1:nrow(BLCA_Subty)){
  BLCA_Subty$`Days until death`[i]<-max(BLCA_Subty$`Days until death`[i],BLCA_Subty$`Days to last followup`[i])
}
index = which(BLCA_Subty$`Days until death`!="0")
BLCA_Subty = BLCA_Subty[index,]
proportions =proportions[index,]
BLCA_Subty$`Days until death` = as.numeric(BLCA_Subty$`Days until death`)
allData = BLCA_Subty            #402 samples
allData$`Vital status`= allData$`Vital status`== "Dead"

## draw the survival chart
for(j in c(6,9) ){
  thisName = colnames(proportions)[j]
  allData$group = "Low"
  #allData$group[proportions[,thisName] > median(proportions[,thisName])] = "High"
  quantiles = quantile(proportions[,j], prob = seq(0, 1, length = 6))
  m1 = quantiles[2]
  m2 = quantiles[5]
  allData$group[proportions[,thisName] <= m1] = "group1"
  allData$group[proportions[,thisName] >= m2] = "group2"
  allData01 = subset(allData,group == "group1")
  allData02 = subset(allData,group == "group2") 
  allData0 = rbind(allData01,allData02)      #162 samples
  SurvObj = with(allData0, Surv(`Days until death`, `Vital status`))
  fit = coxph(SurvObj ~ as.factor(group), allData0)
  pval = round(summary(fit)$coefficients[5],4)
  plot(survfit(SurvObj ~ as.factor(group), data = allData0), col=c("red", "blue"),
         mark.time=TRUE, xlab="Days", ylab="Survivial rate", main = paste0("BLCA:",thisName))
  legend("topright", legend=c("Top 20%", "Bottom 20%"),
         col=c("red", "blue"), lty=c(1,1), pch = c("+","+"),bty="n",cex=1.1)
  legend("bottomleft",legend=paste0("p=",pval), bty="n",cex=1.1)
}
dev.off()



