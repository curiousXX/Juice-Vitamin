setwd('YOUR_PATH')
library(readxl)

#----------------------Calculation indicators--------
Dat <- read_excel("ME213_43.xlsx")
Dat = as.data.frame(Dat)
row.names(Dat) = Dat$ID
Dat = Dat[,-1]
Dat$`Health.NLR` <-Dat$Health.中性粒细胞总数.G.L. / Dat$Health.淋巴细胞总数.G.L.
Dat$`Health.SII` <-Dat$Health.中性粒细胞总数.G.L. * Dat$Health.血小板计数值.G.L./ Dat$Health.淋巴细胞总数.G.L.
saveRDS(Dat,'ME213_45.rds')

#-----------------INPUT---------
Dat = readRDS('ME213_45.rds')
Info = read_excel("MEINFO213.xlsx")
Info = as.data.frame(Info)
saveRDS(Info,'MEINFO213.rds')
#-------------------PCA-------------------------
library(ggplot2)
library(factoextra)
site <- which(Info$group2=='AOS')
pca <- prcomp(Dat[site,],scale. = T)
fviz_eig(pca, addlabels = TRUE)
value <- as.factor(as.character(Info$phase)[site])
fviz_pca_ind(pca, col.ind=value, mean.point=T, label='none',
             addEllipses = T,
             legend.title="Phase")

pca <- prcomp(Dat[-site,],scale. = T)
fviz_eig(pca, addlabels = TRUE)
value <- as.factor(as.character(Info$phase)[-site])
fviz_pca_ind(pca, col.ind=value, mean.point=T, label='none',
             addEllipses = T,
             legend.title="Phase")

site <- which(Info$group1=='Juice')
pca <- prcomp(Dat[site,],scale. = T)
fviz_eig(pca, addlabels = TRUE)
value <- as.factor(as.character(Info$phase)[site])
fviz_pca_ind(pca, col.ind=value, mean.point=T, label='none',
             addEllipses = T,
             legend.title="Phase")

site <- which(Info$group1=='VE')
subdata <- kill_0(Dat[site,])
pca <- prcomp(subdata,scale. = T)
fviz_eig(pca, addlabels = TRUE)
value <- as.character(Info$phase)[fsite(row.names(subdata),Info,'ID')]
fviz_pca_ind(pca, col.ind=value, mean.point=T, label='none',
             addEllipses = T,
             legend.title="Phase")

site <- which(Info$group1=='GSE')
pca <- prcomp(Dat[site,],scale. = T)
fviz_eig(pca, addlabels = TRUE)
value <- as.factor(as.character(Info$phase)[site])
fviz_pca_ind(pca, col.ind=value, mean.point=T, label='none',
             addEllipses = T,
             legend.title="Phase")



#----------------------adonis2---------
library(vegan)
s1 <- which(Info$phase2==0)
anvo=adonis2(formula = Dat[s1,]~group1,data = Info[s1,])
anvo
s1 <- which(Info$phase2==1)
adonis2(formula = Dat[s1,]~group2,data = Info[s1,])

s1 <- which(Info$group2=='AOS')
adonis2(formula = Dat[s1,]~phase,data = Info[s1,])

s1 <- which(Info$group2=='Placebo')
adonis2(formula = Dat[s1,]~phase,data = Info[s1,])

s1 <- which(Info$group1=='Juice')
adonis2(formula = Dat[s1,]~phase,data = Info[s1,])

s1 <- which(Info$group1=='GSE')
adonis2(formula = Dat[s1,]~phase,data = Info[s1,])

s1 <- which(Info$group1=='VE')
adonis2(formula = Dat[s1,]~phase,data = Info[s1,])

subdata <- Dat[-which(Info$group2=='Placebo'),]
subinfo <- Info[-which(Info$group2=='Placebo'),]
s1 <- which(subinfo$phase2==0)
adonis2(formula = subdata[s1,]~group1,data = subinfo[s1,])

s1 <- which(subinfo$phase2==1)
adonis2(formula = subdata[s1,]~group1,data = subinfo[s1,])




#------------------------------------test 0 vs 1------------------------------------------
source('F_know_data.R')
g1 <- c('Juice','VE','GSE','Placebo')
g2 <- c('AOS','Placebo')

site <- which(Info$phase!='Withdrawal')
subdata <- Dat[site,];subinfo <- Info[site,]
fea <- names(table(subinfo$Name)[table(subinfo$Name)==2])
site <- fsite(fea,subinfo,'Name')
subdata <- subdata[site,];subinfo <- subinfo[site,]

result <- data.frame(group = 1:(4*ncol(subdata)),feature=NA,p_value=NA,log2FC=NA)
k <- 1
for(i in 1:length(g1)){
  site <- which(subinfo$group1==g1[i])
  sub2info <- subinfo[site,]
  sub2data <- subdata[site,]
  
  for(j in 1:ncol(subdata)){
    v1 <- sub2data[sub2info$phase2==0,j]
    v2 <- sub2data[sub2info$phase2==1,j]
    temp <- wilcox.test(v1,v2,paired =T)
    result$group[k] <- g1[i]
    result$feature[k] <-colnames(sub2data)[j] 
    result$p_value[k] <- temp$p.value
    result$log2FC[k] <- log2(sum(v2)/sum(v1))
    k <- k+1
  }
}

site <- which(result$p_value=='NaN')
result <- result[-site,]
result$q_value <- NA
for(i in 1:length(g1)){
  site <- which(result$group==g1[i])
  value <- p.adjust(result$p_value[site],method = 'fdr')
  result$q_value[site] <-value 
}
saveRDS(result,'ME_test_result.rds')
df <- result[result$p_value<0.1,]
site <- which(grepl('Health',result$feature))
subresult <- result[site,]

for(i in 1:length(g1)){
  site <- which(subresult$group==g1[i])
  value <- p.adjust(subresult$p_value[site],method = 'fdr')
  subresult$q_value[site] <-value 
}
saveRDS(subresult,'ME_test_Health_only.rds')

library(ggplot2)
Info$value <- Dat$Health.低密度脂蛋白.mmol.L.
sum(Info$ID==row.names(Dat))
subinfo <- Info[Info$phase2!=2,]
df <- data.frame(table(subinfo$name2))
fea <- df$Var1[df$Freq==1]
subinfo <- subinfo[-fsite(fea,subinfo,'name2'),]
pdf('LDL all timepoint.pdf',width = 6,height = 4)
ggplot(subinfo,aes(group1,value,fill=phase))+
  geom_boxplot(outlier.color = 'grey')+
  scale_fill_manual(values= c('#4DBBD5','#8491B4'))+
  theme_classic()+ylab('LDL-C mmol/L')+xlab('')
dev.off()
#---------------------------test 1 vs 2----------------------
site <- which(Info$phase!='Baseline')
subdata <- Dat[site,];subinfo <- Info[site,]
fea <- names(table(subinfo$Name)[table(subinfo$Name)==2])
site <- fsite(fea,subinfo,'Name')
subdata <- subdata[site,];subinfo <- subinfo[site,]

result <- data.frame(group = 1:(4*ncol(subdata)),feature=NA,p_value=NA,log2FC=NA)
k <- 1
for(i in 1:length(g1)){
  site <- which(subinfo$group1==g1[i])
  sub2info <- subinfo[site,]
  sub2data <- subdata[site,]
  
  for(j in 1:ncol(subdata)){
    v1 <- sub2data[sub2info$phase2==1,j]
    v2 <- sub2data[sub2info$phase2==2,j]
    temp <- wilcox.test(v1,v2,paired =T)
    result$group[k] <- g1[i]
    result$feature[k] <-colnames(sub2data)[j] 
    result$p_value[k] <- temp$p.value
    result$log2FC[k] <- log2(sum(v2)/sum(v1))
    k <- k+1
  }
}

site <- which(result$p_value=='NaN')
if(length(site)>0){result <- result[-site,]}
result$q_value <- NA
for(i in 1:length(g1)){
  site <- which(result$group==g1[i])
  value <- p.adjust(result$p_value[site],method = 'fdr')
  result$q_value[site] <-value 
}
saveRDS(result,'ME_test_result_1_2.rds')
df <- result[result$p_value<0.1,]
site <- which(grepl('Health',result$feature))
subresult <- result[site,]

for(i in 1:length(g1)){
  site <- which(subresult$group==g1[i])
  value <- p.adjust(subresult$p_value[site],method = 'fdr')
  subresult$q_value[site] <-value 
}
saveRDS(subresult,'ME_test_Health_only_1_2.rds')

#---------------------------test 0 vs 2----------------------
site <- which(Info$phase!='Intervention')
subdata <- Dat[site,];subinfo <- Info[site,]
fea <- names(table(subinfo$Name)[table(subinfo$Name)==2])
site <- fsite(fea,subinfo,'Name')
subdata <- subdata[site,];subinfo <- subinfo[site,]

result <- data.frame(group = 1:(4*ncol(subdata)),feature=NA,p_value=NA,log2FC=NA)
k <- 1
for(i in 1:length(g1)){
  site <- which(subinfo$group1==g1[i])
  sub2info <- subinfo[site,]
  sub2data <- subdata[site,]
  
  for(j in 1:ncol(subdata)){
    v1 <- sub2data[sub2info$phase2==0,j]
    v2 <- sub2data[sub2info$phase2==2,j]
    temp <- wilcox.test(v1,v2,paired =T)
    result$group[k] <- g1[i]
    result$feature[k] <-colnames(sub2data)[j] 
    result$p_value[k] <- temp$p.value
    result$log2FC[k] <- log2(sum(v2)/sum(v1))
    k <- k+1
  }
}

site <- which(result$p_value=='NaN')
if(length(site)>0){result <- result[-site,]}
result$q_value <- NA
for(i in 1:length(g1)){
  site <- which(result$group==g1[i])
  value <- p.adjust(result$p_value[site],method = 'fdr')
  result$q_value[site] <-value 
}
saveRDS(result,'ME_test_result_0_2.rds')
df <- result[result$p_value<0.1,]
site <- which(grepl('Health',result$feature))
subresult <- result[site,]

for(i in 1:length(g1)){
  site <- which(subresult$group==g1[i])
  value <- p.adjust(subresult$p_value[site],method = 'fdr')
  subresult$q_value[site] <-value 
}
saveRDS(subresult,'ME_test_Health_only_0_2.rds')
#-------------------------------RTL-----------------------------------
Dat <- readRDS('ME213_126.rds')
df <- readRDS('ME_test_result.rds')
df$comp <- 'B_I'
result <- df[df$feature=='Vitamin.维生素E检测值.ng.mL.',]

df <- readRDS('ME_test_result_1_2.rds')
df$comp <- 'I_W'
temp <- df[df$feature=='Vitamin.维生素E检测值.ng.mL.',]
result <- rbind(result,temp)

df <- readRDS('ME_test_result_0_2.rds')
df$comp <- 'B_W'
temp <- df[df$feature=='Vitamin.维生素E检测值.ng.mL.',]
result <- rbind(result,temp)


df <- readRDS('ME_test_Health_only.rds')
subdf <- df[df$q_value<0.05,]
temp <- data.frame(table(subdf$group,subdf$feature))
df2 <- reshape2::dcast(temp,Var1~Var2)
df2 <- t(df2);colnames(df2) <- df2[1,];df2 <- df2[-1,]
df2[,1] <- as.numeric(df2[,1]);df2[,2] <- as.numeric(df2[,2]);df2[,3] <- as.numeric(df2[,3])
df2$sum <- rowSums(df2)
fea <- row.names(df2)

result<- readRDS('ME_test_Health_only.rds')
result <- result[which(result$feature %in% fea),]
result <- result[,-3];colnames(result)[3:4] <- paste0(colnames(result)[3:4],'_1_2')

temp <- readRDS('ME_test_Health_only_0_2.rds')
temp <- temp[which(temp$feature %in% fea),]
for(i in unique(temp$group)){
  site <- which(temp$group==i)
  temp$q_value[site] <- p.adjust(temp$p_value[site])
}
temp <- temp[,-3];colnames(temp)[3:4] <- paste0(colnames(temp)[3:4],'_0_2')
sum(result$feature==temp$feature);sum(result$group==temp$group)
result <- cbind(result,temp[,3:4])

temp <- readRDS('ME_test_result_1_2.rds')
site <- which(grepl('Health',temp$feature))
temp <- temp[site,]
length(unique(temp$feature))
temp <- temp[which(temp$feature %in% fea),]
for(i in unique(temp$group)){
  site <- which(temp$group==i)
  temp$q_value[site] <- p.adjust(temp$p_value[site])
}


temp <- temp[,-3];colnames(temp)[3:4] <- paste0(colnames(temp)[3:4],'_1_2')
sum(result$feature==temp$feature);sum(result$group==temp$group)
result <- cbind(result,temp[,3:4])
write.csv(result,'Table S3.csv',quote = F,row.names = F)
library(ggplot2)
Info$value <- Dat$Telomere_length
sum(Info$ID==row.names(Dat))
subinfo <- Info[Info$phase2!=2,]
df <- data.frame(table(subinfo$name2))
fea <- df$Var1[df$Freq==1]
subinfo <- subinfo[-fsite(fea,subinfo,'name2'),]
pdf('relative telomere length all timepoint2.pdf',width = 6,height = 4)
ggplot(subinfo,aes(group1,value,fill=phase))+
  geom_boxplot(outlier.color = 'grey')+
  scale_fill_manual(values= c('#4DBBD5','#8491B4'))+
  theme_classic()+ylab('relative telomere length')+xlab('')#+ylim(c(0,1.0))
dev.off()

subinfo <- subinfo[order(subinfo$name2),]
t1 <- subinfo[subinfo$phase2==0,]
t2 <- subinfo[subinfo$phase2==1,]
t1$value2 <- t2$value-t1$value
sub2info <- t1
v1 <- sub2info$value2[sub2info$group1=='Juice']
v2 <- sub2info$value2[sub2info$group1=='VE']
v3 <- sub2info$value2[sub2info$group1=='GSE']
v4 <- sub2info$value2[sub2info$group1=='Placebo']

wilcox.test(v2,v1)
pdf('relative telomere length detla.pdf',width = 4,height = 8)
ggplot(sub2info,aes(group1,value2,fill=group1))+
  geom_boxplot(outlier.color = 'grey')+
  geom_hline(aes(yintercept = 0),linetype = 'dotted',color='black')+
  scale_fill_manual(values= c('#AD002A','#42B540','#925E9F','#ADB6B6'),guide='none')+
  theme_classic()+ylab('relative telomere length')+xlab('')+ylim(c(-0.5,0.8))
dev.off()



df <- data.frame(table(Info$name2))
fea <- df$Var1[df$Freq!=3]
subinfo <- Info[-fsite(fea,Info,'name2'),]
subdata <- Dat[subinfo$ID,]
df<- readRDS('ME_test_Health_only.rds')
fea <- unique(df$feature[df$q_value<0.05])
fea <- fea[order(fea)]
fea <- fea[c(4,1:3,5:19)]

n <- ceiling(sqrt(length(fea)))
site <- which(subinfo$group1!='Placebo')
subinfo <- subinfo[site,];subdata <- subdata[site,]
fea_names <- c('LDL-C mmol/L','GGT U/L','A/G ratio','Albumin g/L',
               'RBC 10*12/L','Hematocrit L/L','IBIL μmol/L','MCH pg',
               'MCHC g/L','MPV fL','Globulin g/L','Eosinophils %',
               'AST U/L','AST/ALT ratio','Creatinine μmol/L','PCT %',
               'DBIL μmol/L','TC mmol/L','TP g/L')
site <- order(fea_names)
fea <- fea[site];fea_names <- fea_names[site]
site <- c(12,1:11,13:19)
fea <- fea[site];fea_names <- fea_names[site]

n1=4;n2=5
pdf('Figure S1.pdf',width = 4*n1,height = 3*n2)
grid::grid.newpage() 
grid::pushViewport(grid::viewport(layout = grid::grid.layout(n2,n1)))
s1 <- s2 <- 1
for(i in 1:length(fea)){
  subinfo$value <- subdata[,fea[i]]
  p <- ggplot(subinfo)+
    geom_boxplot(aes(phase,value,fill=phase),outlier.color = 'white',show.legend =F,size=0.2)+
    scale_fill_manual(values= c('#4DBBD5','#8491B4','#3C5488'))+
    facet_wrap(~group1,ncol =4)+
    geom_point(aes(x=phase,y=value),color='grey',size=0.5)+
    geom_line(aes(x=phase,y=value,group=name2),color='grey',size=0.2)+
    theme_bw()+ylab(fea_names[i])+xlab('')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  print(p,vp = vplayout(s1,s2))
  if(s2==n1){
    s1 <- s1+1
    s2=1
  }
  else{s2=s2+1}                       
  
}
dev.off()

sum(row.names(Dat)==Info$ID)
Info$value <- Dat$Vitamin.维生素E检测值.ng.mL.
subinfo <- Info[Info$group1=='VE',]
temp <- data.frame(table(subinfo$name2))

pdf('SVE_boxplot.pdf',width = 3,height = 5)
ggplot(subinfo)+
  geom_boxplot(aes(phase,value,fill=phase),outlier.color = 'grey',show.legend =F)+
  scale_fill_manual(values= c('#4DBBD5','#8491B4','#3C5488'))+
  theme_classic()+ylab('serum vitamin E ng/mL')+xlab('')
dev.off()

sub2info <- subinfo[subinfo$phase2!=2,]
temp <- data.frame(table(sub2info$name2))
sub2info <- sub2info[fsite(temp$Var1[temp$Freq==2],sub2info,'name2'),]
sub2info$LDL <- Dat[sub2info$ID,'Health.低密度脂蛋白.mmol.L.']
df <- sub2info[,c('name2','phase','value')]
df <- dcast(df,name2~phase)
colnames(df) <- c('name','B_VE','I_VE')
result <- df
df <- sub2info[,c('name2','phase','LDL')]
df <- dcast(df,name2~phase)
colnames(df) <- c('name','B_LDL','I_LDL')
result <- merge(result,df)
result$LDL <- (result$I_LDL-result$B_LDL)/result$B_LDL*100
result$VE <- (result$I_VE-result$B_VE)/result$B_VE*100

pdf('serum vitamin E.pdf',width = 5,height = 5)
ggplot()+  
  geom_point(data = sub2info,mapping = aes(value,LDL,color=phase),size=3,show.legend = T)+
  scale_color_manual(values= c('#4DBBD5','#8491B4'))+
  geom_segment(data=result, mapping = aes(x=B_VE, y=B_LDL, xend=I_VE, yend=I_LDL), color='black',arrow = arrow(length=unit(0.2, "cm"))) + 
  theme_classic()+xlab('serum vitamin E ng/mL')+ylab('LDL-C mmol/L')
dev.off()

#-----------------lipid----------------
library(ggplot2)
Dat <- readRDS('ME213_126.rds')
Info <- readRDS('MEINFO213.rds')

fea <- c('Health.低密度脂蛋白.mmol.L.','Health.高密度脂蛋白.mmol.L.',
         'Health.甘油三酯值.mmol.L.','Health.总胆固醇值.mmol.L.')
fea_name <- c('LDL-C','HDL-C','TG','TC')

for(i in 1:length(fea)){
  Info$value <- Dat[,fea[i]]
  p <- ggplot(Info, aes(x=group1, y=value, fill=phase)) + 
    scale_fill_manual(values= c('#B09C85','#00468B','#ED0000'))+
    geom_boxplot(outlier.colour='grey')+theme_classic()+
    xlab('')+ylab(paste0(fea_name[i],' mmol/L'))
  pdf(paste0(fea_name[i],'_boxplot.pdf'),width = 6,height = 7)
  print(p)
  dev.off()
}

#-----------------------Metabolic changes before and after intervention-------------------
my_group = as.character(unique(Info$group1))
value = stringr::str_remove_all(colnames(Dat),pattern = '\\.')
colnames(Dat) = value

Dat = Dat[,-62]
df = data.frame(fea = colnames(Dat))
my_fea = colnames(Dat)[c(1:22,25:42,124:125)]
temp = know_data(Dat[my_fea],draw = T)
result = data.frame(matrix(NA,nrow = 42*4,ncol = 4))
colnames(result) = c('Group','Feature','Baseline Value','Intervention Value')
fea = c('Health.γ.谷氨酰转肽酶检测值.U.L.','Health.血清谷草转氨酶值.U.L.','Health.血清谷丙转氨酶值.U.L.',
        'Health.甘油三酯值.mmol.L.',
        'Health.直接胆红素.μmol.L.','Health.血清总胆红素检测值.μmol.L.','Health.间接胆红素.μmol.L.')
k=1
for(i in 1:length(my_group)){
  result$Group[k:(k+41)] <- my_group[i]
  result$Feature[k:(k+41)] <- my_fea
  for(j in 1:42){
    v1 = Dat[Info$group1==my_group[i] & Info$phase=='Baseline',my_fea[j]]
    v2 = Dat[Info$group1==my_group[i] & Info$phase=='Intervention',my_fea[j]]
    result$`Baseline Value`[k] <- paste0(round(mean(v1),1),'±',round(sd(v1),1))
    result$`Intervention Value`[k] <- paste0(round(mean(v2),1),'±',round(sd(v2),1))
    if(my_fea[j] %in% fea){
      result$`Baseline Value`[k] <-paste0(round(median(v1),1),' (',round(quantile(v1)[2],1),',',round(quantile(v1)[4],1),')')
      result$`Intervention Value`[k] <-paste0(round(median(v2),1),' (',round(quantile(v2)[2],1),',',round(quantile(v2)[4],1),')')
    }
    k=k+1
  }
}

exdata = readRDS('ME_test_Health_only.rds')
exdata = exdata[,-3]
row.names(result) = paste0(result$Group,'_',result$Feature)
row.names(exdata) = paste0(exdata$group,'_',exdata$feature)
length(intersect(result$Feature,exdata$feature))
exdata = exdata[row.names(result),]

sum(row.names(result)==row.names(exdata))
result = cbind(result,exdata[,3:4])
result$log2FC = round(result$log2FC,3)
result$q_value = round(result$q_value,3)

openxlsx::write.xlsx(result,'Metabolic changes.xlsx')
