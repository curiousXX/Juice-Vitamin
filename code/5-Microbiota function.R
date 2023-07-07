setwd('YOUR_PATH')
source('F_know_data.R')

df <- read.csv("VJ_meta_wilcox_st.csv")

fea <- unique(df$fea)
fea1 <- unique(df$fea[df$F1>0])
fea2 <- unique(df$fea[df$F1<0])
exdata <- readRDS('ko_list.rds')
micro_anno <- readRDS('micro_anno_v2.rds')
fea1 <- c(fea1,'GUT_GENOME143094','GUT_GENOME001779')

exdata <- exdata[exdata$bact_name %in% fea1,]



library(readr)
library(readxl)
anno <- read_excel("lipid_related_ko.xlsx")

exdata <- exdata[exdata$KO %in% anno$KO,]
micro_anno <- micro_anno[ micro_anno$Genome_accession %in% fea1,]
micro_anno$species[micro_anno$Genome_accession=='GUT_GENOME014039'] <- 'Faecalibacterium sp014039'
micro_anno$species[micro_anno$Genome_accession=='GUT_GENOME154740'] <- 'Faecalibacterium sp154740'
exdata$bacteria <- micro_anno$species[fsite(exdata$bact_name,micro_anno,'Genome_accession')]

exdata$value <- 1
exdata$value[fsite(fea2,exdata,'bact_name')] <- -1
Dat <- merge(exdata,anno, all = T)

Dat <- Dat[,-c(2,3)]
Dat <- unique(Dat)
Dat$annotation <- paste0(Dat$KO,': ',Dat$annotation)
saveRDS(Dat,'lipid_ko_in_strains.rds')
library(pheatmap)
library(reshape2)
fea <- unique(Dat$class)
for(i in fea[c(1:3,5,9)]){
  subdata <- Dat[Dat$class==i,c(3,4,5)]
  subdata <- unique(subdata)#cbw
  if(nrow(subdata)<3){
    print(paste0(i,'less than 3'))
    next()}
  df <- dcast(subdata,bacteria~annotation,fill = 0)
  df <- na.omit(df)
  row.names(df) <- df[,1];df <- df[,-1]
  df <- kill_0(df)
  pdf(paste0('pheatmap_lipid_strains_ko_',i,'.pdf'),width = 15,height = 15)
  pheatmap(t(df),cellwidth = 10,cellheight = 7,fontsize = 6,
           color = c('#4DBBD5','white','#E64B35'),legend = F,angle_col = 45)
  dev.off()
}


for(i in fea[c(1:3,5,9)]){
  subdata <- Dat[Dat$class==i,c(3,4,5)]
  subdata <- unique(subdata)
  if(nrow(subdata)<3){
    print(paste0(i,'less than 3'))
    next()}
  df <- dcast(subdata,bacteria~annotation,fill = 0)
  df <- na.omit(df)
  row.names(df) <- df[,1];df <- df[,-1]
  df <- kill_0(df)
  pd=T
  if(ncol(df)<3){pd <- F}
  pdf(paste0('pheatmap_lipid_strains_ko_',i,' 2.pdf'),width = 15,height = 15)
  pheatmap(t(df),cellwidth = 10,cellheight = 7,fontsize = 6,cluster_rows = pd,
           color = c('white','#E64B35'),legend = F,angle_col = 45)
  dev.off()
}

subdata <- Dat[Dat$class=='cholesterol metabolism',c(3,4,5)]
subdata <- unique(subdata)
df <- dcast(subdata,bacteria~annotation,fill = 0)
df <- na.omit(df)
row.names(df) <- df[,1];df <- df[,-1]
pdf('pheatmap_lipid_strains_ko_TC.pdf',width = 15,height = 15)
pheatmap(t(df),cellwidth = 10,cellheight = 7,fontsize = 6,
         color = c('#4DBBD5','white','#E64B35'),legend = F,angle_col = 45)
dev.off()

subdata <- Dat[Dat$class=='Lys-type peptidoglycan biosynthesis',c(3,4,5)]
subdata <- unique(subdata)
df <- dcast(subdata,bacteria~annotation,fill = 0) #cwb
df <- na.omit(df)
row.names(df) <- df[,1];df <- df[,-1]
pdf('pheatmap_lipid_strains_ko_PG.pdf',width = 15,height = 15)
pheatmap(t(df),cellwidth = 10,cellheight = 7,fontsize = 6,
         color = c('white','#E64B35'),legend = F,angle_col = 45)
dev.off()

subdata <- Dat[Dat$class=='Exopolysaccharide biosynthesis',c(3,4,5)]
subdata <- unique(subdata)
df <- dcast(subdata,bacteria~annotation,fill = 0)
df <- na.omit(df)
row.names(df) <- df[,1];df <- df[,-1]
pdf('pheatmap_lipid_strains_ko_EPS.pdf',width = 15,height = 15)
pheatmap(t(df),cellwidth = 10,cellheight = 7,fontsize = 6,
         color = c('white','#E64B35'),legend = F,angle_col = 45)
dev.off()


subdata <- Dat[Dat$class=='Indole production',c(3,4,5)]
subdata <- unique(subdata)

df <- dcast(subdata,bacteria~annotation,fill = 0)
df <- na.omit(df)
row.names(df) <- df[,1];df <- df[,-1]
pdf('pheatmap_lipid_strains_ko_SCFA.pdf',width = 15,height = 15)
pheatmap(t(df),cellwidth = 10,cellheight = 7,fontsize = 6,cutree_cols = 2,
         color = c('#4DBBD5','white','#E64B35'),legend = F,angle_col = 45)
dev.off()







