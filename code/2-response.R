source('F_know_data.R')
setwd('YOUR_PATH')
library(readxl)
library(ARTool)
topx <- function(value,x=3){
  y <- order(value)
  s1 <- y[1:3]
  y <- order(value,decreasing = T)
  s2 <- y[1:3]
  return(c(s1,s2))
}


#---------------------------------data process---------------------
Dat <- readRDS('ME213_45.rds')
Info <- readRDS('MEINFO213.rds')
t1 <- read_excel("Info112.xlsx")
t1 = as.data.frame(t1)
saveRDS(t1,'Info112.rds')

fea <- intersect(Info$id2,t1$id2)
site <- t1$seq2[fsite(fea,t1,'id2')]
Info <- Info[fsite(fea,Info,'id2'),]
site <- fsite(c('Juice','VE'),Info,'group1')
Info <- Info[site,]
colnames(Info)[4] <- 'name'
intersect(colnames(Info),colnames(t1))
Info <- merge(Info,t1[,-c(2,4,17,18,7)],all.x = T,by='id2')
Dat <- Dat[Info$ID,]
saveRDS(Info,'MEINFO47.rds')
saveRDS(Dat,'ME47_45.rds')

Dat <- read_excel("D112_1098.xlsx")
Dat = as.data.frame(Dat)
row.names(Dat) <- Dat$ID
Dat = Dat[,-1]
saveRDS(Dat,'D112_1098.rds')
#---------------------------------VE--------------------------
Info <- readRDS('MEINFO47.rds')
Dat <- readRDS('ME47_45.rds')

Info$value <- Dat$Health.低密度脂蛋白.mmol.L.
site <- which(Info$group1=='VE' & Info$phase!='Withdrawal' & Info$id2!='0244')

pdf('VE_with_GUTDATA_LDL.pdf',width = 3,height = 4)
ggplot(Info[site,])+
  geom_boxplot(aes(x=phase,y=value,color=phase))+ guides(color='none')+
  scale_color_manual(values= c('#4DBBD5','#8491B4'))+
  geom_point(aes(x=phase,y=value),color='grey')+
  geom_line(aes(x=phase,y=value,group = name3),color='grey')+
  theme_classic()+ylab('LDL-C')+xlab('')

dev.off()

df <- Info[site,]
fea <- unique(df$name3)
result <- data.frame(names=fea,change=NA,change_ratio=NA)
for(i in 1:length(fea)){
  v1 <- df$value[df$name3==fea[i] & df$phase2==0]
  v2 <- df$value[df$name3==fea[i] & df$phase2==1]
  result$change[i] <- round(v2-v1,3)
  result$change_ratio[i] <- round((v2-v1)/v1,3)
}



site <- topx(result$change_ratio)
write.csv(result,'VE_LDL_change_ratio.csv',quote = F,row.names = F)
result <- result[site,]
result$group <- c(rep('top',3),rep('bottom',3))
temp <- Info[Info$phase2!=2,]
id <- temp$seq2[fsite(result$names,temp,'name3')]

mean(result$change_ratio[1:3]*100);sd(result$change_ratio[1:3]*100)
mean(result$change_ratio[4:6]*100);sd(result$change_ratio[4:6]*100) #cbw

#---------------------------------------------

Dat = readRDS('D112_1098.rds')
Dat <- Dat[id,]
Dat <- kill_0(Dat)
site <- c()
for(i in 1:ncol(Dat)){
  value <- Dat[,i]
  value <- as.numeric(value>0)
  if(sum(value)<2){site <- c(site,i)}
}
Dat <- Dat[,-site]
df <- data.frame(id =row.names(Dat),group =NA,name =NA,timepoint =NA)
df$name <- Info$name3[fsite(df$id,Info,'seq2')]
df$timepoint <- as.factor(Info$phase[fsite(df$id,Info,'seq2')])
df$group <- as.factor(result$group[fsite(df$name,result,1)])
df$name <- as.factor(df$name)

my_result <- data.frame(feature=colnames(Dat),P_group=NA,P_timepoint=NA,P_interact=NA)
for(i in 1:ncol(Dat)){
  df$value <- Dat[,i]
  m <- art(value ~ group*timepoint + Error(name), data=df)
  temp <- anova(m)
  my_result[i,2:4] <- temp$`Pr(>F)`
}
my_result$Q_group <- p.adjust(my_result$P_group,'fdr')
microbe_anno <- readRDS('micro_anno_v2.rds')
my_result$species <- microbe_anno$species[fsite(my_result$feature,microbe_anno,13)]
my_result$anno <- microbe_anno$Lineage[fsite(my_result$feature,microbe_anno,13)]

write.csv(my_result,'VE_twoway_anova.csv',quote = F,row.names = F)

fea <- my_result$feature[my_result$P_group<0.05]
df2 <- data.frame(feature =NA,group=NA,value=1:(length(fea)*4))
k=1
for(i in 1:length(fea)){
  df2$feature[k:(k+3)] <- fea[i]
  df2$group[k:(k+3)] <-c('top_0','top_1','bottom_0','bottom_1')
  site <- which(df$group=='top' & df$timepoint=='Baseline')
  value <- Dat[site,fea[i]]
  df2$value[k] <- mean(value)
  site <- which(df$group=='top' & df$timepoint=='Intervention')
  value <- Dat[site,fea[i]];df2$value[k+1] <- mean(value)
  site <- which(df$group=='bottom' & df$timepoint=='Baseline')
  value <- Dat[site,fea[i]];df2$value[k+2] <- mean(value)
  site <- which(df$group=='bottom' & df$timepoint=='Intervention')
  value <- Dat[site,fea[i]];df2$value[k+3] <- mean(value)
  k=k+4
}
df2$bacteria <- microbe_anno$species[fsite(df2$feature,microbe_anno,13)] 
library(pheatmap)
ma <- reshape2::dcast(df2[,-1],bacteria ~group)
row.names(ma) <- ma[,1];ma <- ma[,-1]
ma <- ma[,c(3,4,1,2)]
min(ma[ma>0])
pdf('VE_anova_heatmap.pdf',width = 10,height = 10)
pheatmap(log(ma+1e-10,10),cluster_cols  = F,cluster_rows = T,
         main='top bottom',
         labels_col=c('Baseline','Intervention','Baseline','Intervention'),
         cellwidth = 39,cellheight = 10,fontsize = 8,angle_col = 0,
         colorRampPalette(colors = c("#3C5488","white","#DC0000"))(200))
dev.off()
pheatmap(log(ma+1e-10,10),cluster_cols  = F,cluster_rows = F,
         cellwidth = 9,cellheight = 9,fontsize = 8,
         colorRampPalette(colors = c("#3C5488","white","#DC0000"))(200))

#----------------------------------------JUICE----------------
site <- which(Info$group1=='Juice' & Info$phase!='Withdrawal' & Info$id2!='0181')
df <- Info[site,]

pdf('Juice_with_GUTDATA_LDL.pdf',width = 3,height = 4)
ggplot(Info[site,])+
  geom_boxplot(aes(x=phase,y=value,color=phase))+ guides(color='none')+
  scale_color_manual(values= c('#4DBBD5','#8491B4'))+
  geom_point(aes(x=phase,y=value),color='grey')+
  geom_line(aes(x=phase,y=value,group = name3),color='grey')+
  theme_classic()+ylab('LDL-C')+xlab('')
dev.off()


df <- Info[site,]
fea <- unique(df$name3)
result <- data.frame(names=fea,change=NA,change_ratio=NA)
for(i in 1:length(fea)){
  v1 <- df$value[df$name3==fea[i] & df$phase2==0]
  v2 <- df$value[df$name3==fea[i] & df$phase2==1]
  result$change[i] <- round(v2-v1,3)
  result$change_ratio[i] <- round((v2-v1)/v1,3)
}

site <- topx(result$change_ratio)
write.csv(result,'Juice_LDL_change_ratio.csv',quote = F,row.names = F)
result <- result[site,]
result$group <- c(rep('top',3),rep('bottom',3))
temp <- Info[Info$phase2!=2,]
id <- temp$seq2[fsite(result$names,temp,'name3')]

Dat <- readRDS('D112_1098.rds') #cbw
Dat <- Dat[id,]
Dat <- kill_0(Dat)
site <- c()
for(i in 1:ncol(Dat)){
  value <- Dat[,i]
  value <- as.numeric(value>0)
  if(sum(value)<2){site <- c(site,i)}
}
Dat <- Dat[,-site]
df <- data.frame(id =row.names(Dat),group =NA,name =NA,timepoint =NA)
df$name <- Info$name3[fsite(df$id,Info,'seq2')]
df$timepoint <- as.factor(Info$phase[fsite(df$id,Info,'seq2')])
df$group <- as.factor(result$group[fsite(df$name,result,1)])
df$name <- as.factor(df$name)

my_result <- data.frame(feature=colnames(Dat),P_group=NA,P_timepoint=NA,P_interact=NA)
for(i in 1:ncol(Dat)){
  df$value <- Dat[,i]
  m <- art(value ~ group*timepoint + Error(name), data=df)
  temp <- anova(m)
  my_result[i,2:4] <- temp$`Pr(>F)`
}
my_result$Q_group <- p.adjust(my_result$P_group,'fdr')
microbe_anno <- readRDS('micro_anno_v2.rds')
my_result$anno <- microbe_anno$Lineage[fsite(my_result$feature,microbe_anno,13)]
my_result$species <- microbe_anno$species[fsite(my_result$feature,microbe_anno,13)]
write.csv(my_result,'Juice_twoway_anova.csv',quote = F,row.names = F)

fea <- my_result$feature[my_result$P_group<0.05]
df2 <- data.frame(feature =NA,group=NA,value=1:(length(fea)*4))
k=1
for(i in 1:length(fea)){
  df2$feature[k:(k+3)] <- fea[i]
  df2$group[k:(k+3)] <-c('top_0','top_1','bottom_0','bottom_1')
  site <- which(df$group=='top' & df$timepoint=='Baseline')
  value <- Dat[site,fea[i]]
  df2$value[k] <- mean(value)
  site <- which(df$group=='top' & df$timepoint=='Intervention')
  value <- Dat[site,fea[i]];df2$value[k+1] <- mean(value)
  site <- which(df$group=='bottom' & df$timepoint=='Baseline')
  value <- Dat[site,fea[i]];df2$value[k+2] <- mean(value)
  site <- which(df$group=='bottom' & df$timepoint=='Intervention')
  value <- Dat[site,fea[i]];df2$value[k+3] <- mean(value)
  k=k+4
}
df2$bacteria <- microbe_anno$species[fsite(df2$feature,microbe_anno,13)]
library(pheatmap)
ma <- reshape2::dcast(df2[,-1],bacteria ~group)
row.names(ma) <- ma[,1];ma <- ma[,-1]
ma <- ma[,c(3,4,1,2)]
min(ma[ma>0])
pdf('Juice_anova_heatmap.pdf',width = 10,height = 10)
pheatmap(log(ma+1e-10,10),cluster_cols  = F,cluster_rows = T,
         main='top bottom',
         labels_col=c('Baseline','Intervention','Baseline','Intervention'),
         cellwidth = 39,cellheight = 30,fontsize = 8,angle_col=0,
         colorRampPalette(colors = c("#3C5488","white","#DC0000"))(200))
dev.off()