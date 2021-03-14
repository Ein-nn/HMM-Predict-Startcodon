#@ ID: 1730416017
#@ Author: Hangyang Zhang
#@ Date: 2019-12-06 19:27:10
#@ Last Modified time: 2019-12-19 13:39:17
#@ Email: 1730416017@stu.suda.edu.cn
# setwd("E:/Script/genomics/step5/")

# ++++++ Visualize the distribution of TATA ++++++
library(ggplot2); library(grid); library(pROC)
tata <- read.table("./relativeLocation.txt", header = T, sep = "\t")
colnames(tata) <- c("chr","tataSite","relativeSite","-log(score)","pValue","chain","geneid","tataseq")
tata <- subset(tata, select = c("chr","tataSite","relativeSite","-log(score)","pValue"))
# whole genome
position <- rep("inner", dim(tata)[1])
for (i in 1:dim(tata)[1]) {
  if (tata[i,3] > 0 && tata[i,3] < 200) {
    position[i] <- "0~200"
  }else if(tata[i,3] < 1000 && tata[i,3] >= 200){
    position[i] <- "200~1K"
  }else if(tata[i,3] < 10000 && tata[i,3] >= 1000){
    position[i] <- "1K~10K"
  }else if(tata[i,3] >= 10000){
    position[i] <- ">10K"
  }
}
position <- factor(position)
ggplot(data = tata, aes(x = tataSite, y = `-log(score)`)) + 
  geom_point(aes(shape = chr, colour = position)) +
  scale_shape_manual(values = 1:8) + 
  labs(title = "All TATA-box") +
  pdf("./tata_distribution.pdf",7,7)
dev.off()

tata <- tata[which(position!="inner"),]
position <- position[which(position!="inner")]
ggplot(data = tata, aes(x = tataSite, y = `-log(score)`)) + 
  geom_point(aes(shape = chr, colour = position)) +
  scale_shape_manual(values = 1:8) + 
  labs(title = "Upstream TATA-box") +
  pdf("./upstream_tata_distribution.pdf",7,7)
dev.off()

# ++++++ part of chomsome ++++++
ta <- tata[na.omit(match(0:100000, tata$startSite)),]
gene <- read.table("./gene_anno.gff", sep = "\t")[,c(1,4,5)]
colnames(gene) <- c("chr", "start", "end")
gen <- gene[na.omit(match(0:100000, gene$start)),]

chr <- levels(gen$chr); li <- list()
for (i in 1:length(chr)) {
  c <- ggplot(data = subset(ta, ta$chr == chr[i]), 
              aes(x = startSite, y = `-logScore`, colour = `p-value`)) + 
    geom_point(alpha = .8) + 
    geom_vline(xintercept = subset(gen, gen$chr == chr[i])$start, col = "red") +
    geom_vline(xintercept = subset(gen, gen$chr == chr[i])$end) +
    annotate("text", x=105000, y=9, label=paste("chr",chr[i]))
  if (i < length(chr)){
    c <- c + labs(x = "")
  }
  li <- c(li, list(c))
}
grid.newpage()
pushViewport(viewport(layout = grid.layout(8,1)))
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
for (i in 1:length(chr)) {
  print(li[i], vp = vplayout(i,1)) 
}
dev.off()

# ++++++ Visualize the relative location  ++++++
relative_loc <- read.table("./relativeLocation.txt", header = T, sep = "\t")
rl <- relative_loc$relativeSite
rl1 <- min(rl); rl2 <- max(rl); ave <- mean(rl); med <- median(rl); 
# rl2 - rl1 = 34969; ave = 39.6; med = -11; mode = 1 ,count = 70;quan = -817  525
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
mode <- getmode(rl); quan <- quantile(rl)
nor <- c(ave, med, mode, quan[c(2,4)])
hist(rl, freq = T, breaks = 10000, xlab = "relative site", main = "distribution plot")
hist(rl, freq = T, breaks = 10000, xlim = c(0,200),
     xlab = "relative site", main = "distribution plot")
abline(v = nor,lty=3,lwd=3,col=1:length(nor))
legend("topright", 
       lty = 3, lwd = 3, col = 1:length(nor),
       legend = c("ave: 39.6", "med: -11", "mode: 1", "quan1: -817", "quan2: 525")) 

# +++++++++++++++ ROC ++++++++++++++++++++++++++++++++++++++++++++++++++
rela_loc <- read.table("./relativeLocation.txt", header = T, sep = "\t")
colnames(rela_loc) <- c("chr","tataSite","relativeSite","-log(score)","pValue","chain","geneid","tataseq")
# colnames(rela_loc) <- c("chr","ts","rs","scoreM","chain","geneid","seq")
ta <- subset(rela_loc,select = c("chr","relativeSite","-log(score)"))
ptfp <- function(x = 100){
  flag <- rep(0, dim(ta)[1])
  for (i in 1:length(flag)) {
    if (ta[i,2] < x && ta[i,2] > 11) {
      flag[i] <- 1
    }
  }
  re <- cbind(flag, ta[,3])
  return(re)
}

for (i in 5:20) {
  re <- ptfp(i*10)
  roc_new <- roc(re[,1], re[,2])
  if (i == 5) {
    plot(roc_new, col = i, print.auc=TRUE,
         print.auc.x=(1-i/20),print.auc.y=i/20)
  }else{
    plot.roc(roc_new, add = T, col = i, print.auc=TRUE,
             print.auc.x=(1-i/20),print.auc.y=i/20)
  }
}
re <- ptfp(150)
r <- roc(re[,1], re[,2], plot=TRUE, print.thres=TRUE, print.auc=TRUE,
         print.auc.x = 0.2, print.auc.y= 0.2,
         main = 'Set 150-0 bp as positive')

# +++++++++++++++ filtered ++++++++++++++++++++++++++++++++++++++++++++++++++
tata <- read.table("./relativeLocation.txt", header = T, sep = "\t")
colnames(tata) <- c("chr","tataSite","relativeSite","-log(score)","pValue","chain","geneid","tataseq")
cutoff <- as.numeric(r$auc)
ta <- tata[which(tata[,4] <= 4.587),]
write.table(ta,file = "TATAfiltered.txt", row.names = F, col.names = F, 
            sep = "\t", quote = F)
