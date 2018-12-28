################
### function ###
################
library("fields")
library("plyr")
library("limma")
library("maptools")
library("apcluster")
library("Heatplus")
# source("cb.R")
# source("C:/Users/Makoto/OneDrive/Documents/永野�??/RNA-seq/Os2015/Os2015-2016_figure/scripts/ng.Colors.R")
source("TCC.R")
source("edgeR.R")
source("GOanalysis_functions.R")
load("ulg.RAP_170118")

##############
### script ###
##############

# set parameters -----------------------------------------------------
fn.description <- "141103-1_GeneDescription_Osa"
fn.attribute <- "170101_SampleAttribute_Adachi.txt"

rpm.date <- "170221"
species <- "Osa_Adachi"

fn.rawcnt <- sprintf("%s-2_rawcnt_%s", rpm.date, species)
fn.rpm <- sprintf("%s-3_rpm_%s", rpm.date, species)

fn.gladachi <- "170510_GenelistAdachi.txt"


# load data ----------------------------------
#load sample attribute 
fn <- sprintf("%s", fn.attribute)
at <- read.delim(fn, header=T, as.is=T)
#rownames(at) <- sprintf("%s_%03d", at[,"library"], at[,"index"])

#load transcript description 
fn <- sprintf("%s", fn.description)
load(fn)

#set data read row 
drr <- des[,"NormalizationGroup"]=="data"
names(drr) <- rownames(des)

# load expression data table 
load(file=sprintf("%s", fn.rawcnt))

#load Phtosynthesis gene list of Dr. Adachi 
fn <- sprintf("%s", fn.gladachi)
gladachi <- read.delim(fn, header=T, as.is=T)

# sample, gene criteria -------------------------------------
plot(sort(colSums(rawcnt)))
use.sample <- colSums(rawcnt)>10^5
at = at[use.sample,]
sum(use.sample) # 95 / 96
rawcnt = rawcnt[,-8]
at = at[-8,]

#Calulate RPM
rpm = t(t(rawcnt)/colSums(rawcnt[des$NormalizationGroup=="data",])*10^6)
log2rpm = log2(rpm+0.1)

#Caluculate ave, sd and se
condition = sort(unique(at$AfterLightInduction))
line = unique(at$LineName)
log2rpm.ave.kos = matrix(0, nrow = nrow(log2rpm), ncol = length(condition))
rownames(log2rpm.ave.kos) = rownames(log2rpm)
colnames(log2rpm.ave.kos) = condition
log2rpm.ave.tak = log2rpm.ave.kos
log2rpm.sd.kos = log2rpm.ave.kos
log2rpm.sd.tak = log2rpm.ave.kos
log2rpm.se.kos = log2rpm.ave.kos
log2rpm.se.tak = log2rpm.ave.kos
for(i in 1:length(condition)){
  tmp1 = at[at$AfterLightInduction==condition[i],]
  tmp1 = tmp1[tmp1$LineName==line[1],]$sampleID
  tmp1 = log2rpm[,tmp1]
  ave = rowMeans(tmp1)
  sd = apply(tmp1,1,FUN = sd)
  log2rpm.ave.kos[,i] = ave
  log2rpm.sd.kos[,i] = sd
  log2rpm.se.kos[,i] = sd/sqrt(ncol(tmp1))
  
  tmp2 = at[at$AfterLightInduction==condition[i],]
  tmp2 = tmp2[tmp2$LineName==line[2],]$sampleID
  tmp2 = log2rpm[,tmp2]
  ave = rowMeans(tmp2)
  sd = apply(tmp2,1,FUN = sd)
  log2rpm.ave.tak[,i] = ave
  log2rpm.sd.tak[,i] = sd
  log2rpm.se.tak[,i] = sd/sqrt(ncol(tmp2))
}
# Preparation of the input for DEG analysis
tmp = unique(sort(at$AfterLightInduction))
condition = NULL
name = NULL
for(i in 1:(length(tmp)-1)){
  condition = c(condition, list(c(0,tmp[i+1])))
  name = c(name,sprintf("%sVS%s",0,tmp[i+1]))
}

# Detection of DEGs along time-course
DEG.result.time.course = c()
for(i in 1:length(condition)){
  cat(sprintf("%s %s\n",condition[[i]][1],condition[[i]][2]))
  tmp1 = at[at$AfterLightInduction==condition[[i]][1],]
  tmp1 = tmp1[tmp1$LineName==line[1],]$sampleID
  tmp1 = rawcnt[,tmp1]
  tmp1 = tmp1[des$NormalizationGroup=="data",]
  tmp2 = at[at$AfterLightInduction==condition[[i]][2],]
  tmp2 = tmp2[tmp2$LineName==line[1],]$sampleID
  tmp2 = rawcnt[,tmp2]
  tmp2 = tmp2[des$NormalizationGroup=="data",]
  data = cbind(tmp1, tmp2)
  param_G1 <- ncol(tmp1)                          #G1群のサンプル数を指???
  param_G2 <- ncol(tmp2)                          #G2群のサンプル数を指???
  filename = sprintf("%s",name[i])
  # edgeR.out = edgeR(data = data, param_G1 = param_G1, param_G2 = param_G2, filename = filename)
  tcc.out.kos = tcc(data = data, param_G1 = param_G1, param_G2 = param_G2, filename = filename)
  
  tmp1 = at[at$AfterLightInduction==condition[[i]][1],]
  tmp1 = tmp1[tmp1$LineName==line[2],]$sampleID
  tmp1 = rawcnt[,tmp1]
  tmp1 = tmp1[des$NormalizationGroup=="data",]
  tmp2 = at[at$AfterLightInduction==condition[[i]][2],]
  tmp2 = tmp2[tmp2$LineName==line[2],]$sampleID
  tmp2 = rawcnt[,tmp2]
  tmp2 = tmp2[des$NormalizationGroup=="data",]
  data = cbind(tmp1, tmp2)
  param_G1 <- ncol(tmp1)                          #G1群のサンプル数を指???
  param_G2 <- ncol(tmp2)                          #G2群のサンプル数を指???
  filename = sprintf("%s",condition[i])
  # edgeR.out = edgeR(data = data, param_G1 = param_G1, param_G2 = param_G2, filename = filename)
  tcc.out.tak = tcc(data = data, param_G1 = param_G1, param_G2 = param_G2, filename = filename)

  tmp = list(kos = tcc.out.kos, tak = tcc.out.tak)
  DEG.result.time.course = c(DEG.result.time.course,list(tmp))
}
names(DEG.result.time.course) = name
save(DEG.result.time.course, file = "fig6_DEGresultTimeCourse")


# Detection of DEGs between kos and tak
DEG.result = c()
for(i in c(0,1,5,10,30,60)){
  tmp1 = at[at$AfterLightInduction==i,]
  tmp1 = tmp1[tmp1$LineName==line[1],]$sampleID
  tmp1 = rawcnt[,tmp1]
  tmp1 = tmp1[des$NormalizationGroup=="data",]
  tmp2 = at[at$AfterLightInduction==i,]
  tmp2 = tmp2[tmp2$LineName==line[2],]$sampleID
  tmp2 = rawcnt[,tmp2]
  tmp2 = tmp2[des$NormalizationGroup=="data",]
  data = cbind(tmp1, tmp2)
  param_G1 <- ncol(tmp1)                          #G1群のサンプル数を指???
  param_G2 <- ncol(tmp2)                          #G2群のサンプル数を指???
  filename = sprintf("%s",condition[i])
  tcc.out = tcc(data = data, param_G1 = param_G1, param_G2 = param_G2, filename = filename)
  tmp = list(tcc.out)
  names(tmp) = i
  DEG.result = c(DEG.result, tmp)
}
DEG.result.list = list()
for(i in 1:length(DEG.result)){
  tmp = DEG.result[[i]]
  list = rownames(tmp)[tmp$q.value<0.05]
  tmp = list(list)
  names(tmp) = names(DEG.result)[i]
  DEG.result.list = c(DEG.result.list,tmp)
}
save(DEG.result,DEG.result.list , file = "kos-tak_DEGresult")

DEG.num = c(sum(DEG.result$`0`$q.value<0.05),sum(DEG.result$`1`$q.value<0.05), sum(DEG.result$`5`$q.value<0.05), sum(DEG.result$`10`$q.value<0.05), sum(DEG.result$`30`$q.value<0.05), sum(DEG.result$`60`$q.value<0.05))
names(DEG.num) = c(0,1,5,10,30,60)
pdf("fig6c.pdf")
plot(names(DEG.num),DEG.num, type = "b", pch = 16, ylim = c(0,3500))
dev.off()

# make summary of DEG analysis
DEG.num.time.course = matrix(0,nrow=2, ncol = length(condition))
colnames(DEG.num.time.course) = name
rownames(DEG.num.time.course) = c("kos","tak")
for(i in 1:length(condition)){
  DEG.num.time.course[1,i] = sum(DEG.result.time.course[[i]]$kos$q.value<0.05)
  DEG.num.time.course[2,i] = sum(DEG.result.time.course[[i]]$tak$q.value<0.05)
}

pdf("fig6b.pdf")
condition1 = c("0VS1","0VS5","0VS10","0VS30","0VS60")
condition2 = c("0VS1","0VS5","0VS10")

plot(c(1,5,10,30,60),DEG.num.time.course[2,condition1],type = "n", xaxt = "n", ylab = "Num. of DEGs")
axis(side = 1, at=c(1,5,10,30,60), labels = condition1)
points(c(1,5,10,30,60),DEG.num.time.course[1,condition1], col = "black", pch = 16)
points(c(1,5,10,30,60),DEG.num.time.course[2,condition1], col = "red", pch = 16)
lines(c(1,5,10,30,60),DEG.num.time.course[1,condition1], col = "black")
lines(c(1,5,10,30,60),DEG.num.time.course[2,condition1], col = "red")

plot(c(1,5,10),DEG.num.time.course[2,condition2],type = "n", xaxt = "n", ylab = "Num. of DEGs")
axis(side = 1, at=c(1,5,10), labels = condition2)
points(c(1,5,10),DEG.num.time.course[1,condition2], col = "black", pch = 16)
points(c(1,5,10),DEG.num.time.course[2,condition2], col = "red", pch = 16)
lines(c(1,5,10),DEG.num.time.course[1,condition2], col = "black")
lines(c(1,5,10),DEG.num.time.course[2,condition2], col = "red")
dev.off()


#Obtaine DEG id
DEG.list.time.course.kos = c()
DEG.list.time.course.tak = c()
for(i in 1:length(DEG.result.time.course)){
  tcc.out.kos = rownames(DEG.result.time.course[[i]]$kos[DEG.result.time.course[[i]]$kos$q.value<0.05,])
  tcc.out.tak = rownames(DEG.result.time.course[[i]]$kos[DEG.result.time.course[[i]]$tak$q.value<0.05,])
  DEG.list.time.course.kos = c(DEG.list.time.course.kos,list(tcc.out.kos))
  DEG.list.time.course.tak = c(DEG.list.time.course.tak,list(tcc.out.tak))
}
names(DEG.list.time.course.kos) = names(DEG.result.time.course)
names(DEG.list.time.course.tak) = names(DEG.result.time.course)

#Caluculate common num of DEGs kos VS tak
DEG.kos = unique(unlist(DEG.list.time.course.kos))
DEG.tak = unique(unlist(DEG.list.time.course.tak))

gnl = unique(c(DEG.kos,DEG.tak))

#pca
pdf("fig6a.pdf", width = 14)
layout(matrix(c(1,2),ncol =2))
pca = prcomp(log2rpm[gnl,], scale = T)
pc.rotation = pca$rotation
#Caluculate ave, sd and se
condition = sort(unique(at$AfterLightInduction))
line = unique(at$LineName)
pc.ave.kos = matrix(0, nrow = nrow(pc.rotation), ncol = length(condition))
colnames(pc.ave.kos) = condition
pc.sd.kos = pc.ave.kos
pc.se.kos = pc.ave.kos
pc.ave.tak = pc.ave.kos
pc.sd.tak = pc.ave.kos
pc.se.tak = pc.ave.kos
for(i in 1:length(condition)){
  tmp1 = at[at$AfterLightInduction==condition[i],]
  tmp1 = tmp1[tmp1$LineName==line[1],]$sampleID
  tmp1 = pc.rotation[tmp1,]
  ave = colMeans(tmp1)
  sd = apply(tmp1,2,FUN = sd)
  pc.ave.kos[,i] = ave
  pc.sd.kos[,i] = sd
  pc.se.kos[,i] = sd/sqrt(ncol(tmp1))
  tmp1 = at[at$AfterLightInduction==condition[i],]
  tmp1 = tmp1[tmp1$LineName==line[2],]$sampleID
  tmp1 = pc.rotation[tmp1,]
  ave = colMeans(tmp1)
  sd = apply(tmp1,2,FUN = sd)
  pc.ave.tak[,i] = ave
  pc.sd.tak[,i] = sd
  pc.se.tak[,i] = sd/sqrt(ncol(tmp1))
}
pc.ave = cbind(pc.ave.kos, pc.ave.tak)
pc.sd = cbind(pc.sd.kos, pc.sd.tak)
pc.se = cbind(pc.se.kos, pc.se.tak)
plot(pc.ave[1,], pc.ave[2,],pch = 16, col = c(rep("black",6), rep("red",6)), type = "n", main = "up+down", xlim = range(c(pc.ave[1,]+pc.sd[1,],pc.ave[1,]-pc.sd[1,])), ylim = range(c(pc.ave[2,]+pc.sd[2,],pc.ave[2,]-pc.sd[2,])))
arrows(pc.ave[1,], pc.ave[2,],pc.ave[1,]+pc.sd[1,],pc.ave[2,],angle=90,length=0.05, col = "gray")
arrows(pc.ave[1,], pc.ave[2,],pc.ave[1,]-pc.sd[1,],pc.ave[2,],angle=90,length=0.05, col = "gray")
arrows(pc.ave[1,], pc.ave[2,],pc.ave[1,],pc.ave[2,]+pc.sd[2,],angle=90,length=0.05, col = "gray")
arrows(pc.ave[1,], pc.ave[2,],pc.ave[1,],pc.ave[2,]-pc.sd[2,],angle=90,length=0.05, col = "gray")
points(pc.ave[1,], pc.ave[2,],pch = 16, col = c(rep("black",6), rep("red",6)))
pointLabel(pc.ave[1,], pc.ave[2,], labels=rep(c("0","1","5","10","30","60")),col = c(rep("black",6), rep("red",6)))

plot(pc.ave[1,], pc.ave[3,],pch = 16, col = c(rep("black",6), rep("red",6)), type = "n", main = "up+down", xlim = range(c(pc.ave[1,]+pc.sd[1,],pc.ave[1,]-pc.sd[1,])), ylim = range(c(pc.ave[3,]+pc.sd[3,],pc.ave[3,]-pc.sd[3,])))
arrows(pc.ave[1,], pc.ave[3,],pc.ave[1,]+pc.sd[1,],pc.ave[3,],angle=90,length=0.05, col = "gray")
arrows(pc.ave[1,], pc.ave[3,],pc.ave[1,]-pc.sd[1,],pc.ave[3,],angle=90,length=0.05, col = "gray")
arrows(pc.ave[1,], pc.ave[3,],pc.ave[1,],pc.ave[3,]+pc.sd[3,],angle=90,length=0.05, col = "gray")
arrows(pc.ave[1,], pc.ave[3,],pc.ave[1,],pc.ave[3,]-pc.sd[3,],angle=90,length=0.05, col = "gray")
points(pc.ave[1,], pc.ave[3,],pch = 16, col = c(rep("black",6), rep("red",6)))
pointLabel(pc.ave[1,], pc.ave[3,], labels=rep(c("0","1","5","10","30","60")),col = c(rep("black",6), rep("red",6)))

dev.off()

#Fig6e
pdf("Fig6e.pdf", height = 5)
OsCKX2 = "Os01g0197700"
range = c(-5,13.5)
  gn = OsCKX2
  cols = rep("#FFFFFF00",6) 
  for(n in 1:6){
    sum = sum(is.element(gn, DEG.result.list[[n]]))
    if(sum>0){
      cols[n] = "black"
    }
  }
  plot(c(0,1,5,10,30,60),as.numeric(log2rpm.ave.kos[gn,]),ylim = range, type = "n", main = gn)
  lines(c(0,1,5,10,30,60),as.numeric(log2rpm.ave.kos[gn,]),lwd = 2)
  arrows(c(0,1,5,10,30,60), log2rpm.ave.kos[gn,],c(0,1,5,10,30,60),log2rpm.ave.kos[gn,]+log2rpm.sd.kos[gn,],angle=90,length=0.05,lwd = 2)
  arrows(c(0,1,5,10,30,60), log2rpm.ave.kos[gn,],c(0,1,5,10,30,60),log2rpm.ave.kos[gn,]-log2rpm.sd.kos[gn,],angle=90,length=0.05,lwd = 2)
  lines(c(0,1,5,10,30,60),as.numeric(log2rpm.ave.tak[gn,], col = "red"),lwd = 2)
  arrows(c(0,1,5,10,30,60), log2rpm.ave.tak[gn,],c(0,1,5,10,30,60),log2rpm.ave.tak[gn,]+log2rpm.sd.tak[gn,],angle=90,length=0.05, col = "red",lwd = 2)
  arrows(c(0,1,5,10,30,60), log2rpm.ave.tak[gn,],c(0,1,5,10,30,60),log2rpm.ave.tak[gn,]-log2rpm.sd.tak[gn,],angle=90,length=0.05, col = "red",lwd = 2)
  points(c(0,1,5,10,30,60), apply(cbind(log2rpm.ave.kos[gn,]+log2rpm.sd.kos[gn,]+1,log2rpm.ave.tak[gn,]+log2rpm.sd.tak[gn,]+1),1,FUN = max), pch = "*",col = cols, cex = 2.5)
dev.off()


## GO analysis
#Parse gene id to obtain go terms of use.gene
use.gene = rownames(rawcnt)[rowSums(rawcnt>0)>0] #33816
ulg2 = ulg[is.element(ulg[,"locus"], use.gene),]

# # fisher's exact test for all GO
# result <- ng.mft(cgt=ulg2, gn.test=ERG) 
# ERG.GO = ng.prepGOtestOutTable(result) #0
# fn <- sprintf("go_analysis/ERG.csv")
# write.csv(ng.prepGOtestOutTable(result), file=fn, quote = FALSE)
# result <- ng.mft(cgt=ulg2, gn.test=LRG)
# LRG.GO = ng.prepGOtestOutTable(result) #22
# fn <- sprintf("go_analysis/LRG.csv")
# write.csv(ng.prepGOtestOutTable(result), file=fn, quote = FALSE)
# 
# 
# ng.GetGOTerms(all.enriched.GO[ERG.LRG.GO.mat[,1] ==0 &ERG.LRG.GO.mat[,2] ==0 &ERG.LRG.GO.mat[,3] ==0 &ERG.LRG.GO.mat[,4] ==1])
# ng.GetGOTerms(all.enriched.GO[ERG.LRG.GO.mat[,1] ==0 &ERG.LRG.GO.mat[,2] ==1 &ERG.LRG.GO.mat[,3] ==0 &ERG.LRG.GO.mat[,4] ==0])
# ng.GetGOTerms(all.enriched.GO[ERG.LRG.GO.mat[,1] ==1 &ERG.LRG.GO.mat[,2] ==1])
# ng.GetGOTerms(all.enriched.GO[ERG.LRG.GO.mat[,3] ==1 &ERG.LRG.GO.mat[,4] ==1])
# 
# GO.enrichment.result = list()
# GO.enrichment = function(x,y){
#   out = list()
#   for(i in 1:length(x)){
#     result <- ng.mft(cgt=ulg2, gn.test=x[[i]])
#     # output results as a csv file
#     fn <- sprintf("go_analysis/%s_%s.csv", names(DEG.result.time.course)[i],y)
#     write.csv(ng.prepGOtestOutTable(result), file=fn, quote = FALSE)
#     tmp = list(ng.prepGOtestOut(result))
#     names(tmp) = sprintf("%s_%s", names(DEG.result.time.course)[i],y)
#     out = c(out, tmp)
#   }
#   return(out)
# }
# 
# GO.kos = GO.enrichment(DEG.list.time.course.kos, "kos")
# # GO.kos.up = GO.enrichment(DEG.list.time.course.kos.up, "kos_up")
# # GO.kos.down = GO.enrichment(DEG.list.time.course.kos.down, "kos_down")
# GO.tak = GO.enrichment(DEG.list.time.course.tak, "tak")
# # GO.tak.up = GO.enrichment(DEG.list.time.course.tak.up, "tak_up")
# # GO.tak.down = GO.enrichment(DEG.list.time.course.tak.down, "tak_down")
# result <- ng.mft(cgt=ulg2, gn.test=kostak.60)
# write.csv(ng.prepGOtestOutTable(result), file="go_analysis/kostak_60.csv", quote = FALSE)
# 
# 
# #Make summary of GO-enrichment
# all.enriched.GO = unique(unlist(c(GO.kos,GO.tak)))
# all.enriched.GO.slim = GO.slim(all.enriched.GO)
# all.enriched.GO.slim.catelorized = GO.category(all.enriched.GO.slim)
# check.go = function(x,all){
#   enriched.GO.mat = matrix(0, nrow = length(all), ncol = 5)
#   colnames(enriched.GO.mat) = as.character(c(1,5,10,30,60))
#   rownames(enriched.GO.mat) = all
#   for(i in 1:length(x)){
#     tmp = x[[i]]
#     if(length(tmp)>0){
#       enriched.GO.mat[tmp,i] = 1
#     }
#   }
#   return(enriched.GO.mat)
# }
# GO.kos.mat = check.go(GO.kos,all.enriched.GO)
# GO.tak.mat = check.go(GO.tak,all.enriched.GO)
# GO.mat = cbind(-GO.kos.mat, GO.tak.mat)
# GO.mat2 = rbind(GO.kos.mat, GO.kos.mat)
# for(i in 1:nrow(GO.mat)){
#   GO.mat2[2*i-1,] = GO.tak.mat[i,]
#   GO.mat2[2*i,] = -GO.kos.mat[i,]
# }
# image.plot(t(GO.mat2))
# #Count GO terms of DEGs at each time.point
# GO.count.time.cource.kos.mat = matrix(0, ncol = 5, nrow = length(all.enriched.GO))
# rownames(GO.count.time.cource.kos.mat) = all.enriched.GO
# colnames(GO.count.time.cource.kos.mat) = c(1,5,10,30,60)
# GO.count.time.cource.tak.mat = GO.count.time.cource.kos.mat
# for(i in 1:ncol(GO.count.time.cource.kos.mat)){
#   GO.count.time.cource.kos.mat[,i] = GO.count(DEG.list.time.course.kos[[i]], all.enriched.GO, ulg)  
# }
# for(i in 1:ncol(GO.count.time.cource.tak.mat)){
#   GO.count.time.cource.tak.mat[,i] = GO.count(DEG.list.time.course.tak[[i]], all.enriched.GO, ulg)  
# }
# 
# 
# write.table(sprintf("%s\t%s",rev(all.enriched.GO.slim.catelorized$BP),ng.GetGOTerms(rev(all.enriched.GO.slim.catelorized$BP))), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n", file = "GO_BP.txt")
# write.table(sprintf("%s\t%s",rev(all.enriched.GO.slim.catelorized$MF),ng.GetGOTerms(rev(all.enriched.GO.slim.catelorized$MF))), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n", file = "GO_MF.txt")
# 
# GO.mat[rev(all.enriched.GO.slim.catelorized$BP),]
# GO.mat[rev(all.enriched.GO.slim.catelorized$MF),]


# cytokinin degradation related genes
GO.cytokinin.deg = "GO:0019139"
gnl.cytokinin.deg = ulg[is.element(ulg[,2], GO.cytokinin.deg),]
pdf("Fig6d.pdf", height = 5)
sum(is.element(c(ERG,LRG),gladachi[,3]))
gnl.KEGG.DEG = c(ERG,LRG)[is.element(c(ERG,LRG),gladachi[,3])]
range = c(-5,13.5)
for(i in 1:nrow(gnl.cytokinin.deg)){
  gn = gnl.cytokinin.deg[i]
  cols = rep("#FFFFFF00",6) 
  for(n in 1:6){
    sum = sum(is.element(gn, DEG.result.list[[n]]))
    if(sum>0){
      cols[n] = "black"
    }
  }
  plot(c(0,1,5,10,30,60),as.numeric(log2rpm.ave.kos[gn,]),ylim = range, type = "n", main = gn)
  lines(c(0,1,5,10,30,60),as.numeric(log2rpm.ave.kos[gn,]),lwd = 2)
  arrows(c(0,1,5,10,30,60), log2rpm.ave.kos[gn,],c(0,1,5,10,30,60),log2rpm.ave.kos[gn,]+log2rpm.sd.kos[gn,],angle=90,length=0.05,lwd = 2)
  arrows(c(0,1,5,10,30,60), log2rpm.ave.kos[gn,],c(0,1,5,10,30,60),log2rpm.ave.kos[gn,]-log2rpm.sd.kos[gn,],angle=90,length=0.05,lwd = 2)
  lines(c(0,1,5,10,30,60),as.numeric(log2rpm.ave.tak[gn,]), col = "red",lwd = 2)
  arrows(c(0,1,5,10,30,60), log2rpm.ave.tak[gn,],c(0,1,5,10,30,60),log2rpm.ave.tak[gn,]+log2rpm.sd.tak[gn,],angle=90,length=0.05, col = "red",lwd = 2)
  arrows(c(0,1,5,10,30,60), log2rpm.ave.tak[gn,],c(0,1,5,10,30,60),log2rpm.ave.tak[gn,]-log2rpm.sd.tak[gn,],angle=90,length=0.05, col = "red",lwd = 2)
  points(c(0,1,5,10,30,60), apply(cbind(log2rpm.ave.kos[gn,]+log2rpm.sd.kos[gn,]+1,log2rpm.ave.tak[gn,]+log2rpm.sd.tak[gn,]+1),1,FUN = max), pch = "*",col = cols, cex = 2.5)
}
dev.off()
