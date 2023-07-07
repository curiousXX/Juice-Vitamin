setwd('YOUR_PATH')
source('F_know_data.R')

Info <- readRDS('Info112.rds')
Info <- Info[which(Info$group2_name=='VE' | Info$group2_name=='Juice'),]
table(Info$group2_name,Info$phase2)

Dat <- readRDS('D112_1098.rds')
Dat <- Dat[Info$seq2,]

#---------------------pcoa---------------------------------------
library(ggplot2)
site <- which(Info$group2_name=='VE')
temp <- my_pcoa(Dat[site,])
df <- temp$df
temp$anno
sum(df$names==Info$seq2[site])
df$group <- Info$phase2[site]
ggplot(df,aes(PCoA1,PCoA2,color=group))+geom_point()+theme_classic()

adonis2(Dat[site,]~phase2,data = Info[site,])

temp <- my_pcoa(Dat)
df <- temp$df
temp$anno
sum(df$names==Info$seq2)
df$group <- Info$phase2
df$group2 <- Info$group2_name
pdf('VJ_PCoA.pdf',width = 6,height = 5)
ggplot(df,aes(PCoA1,PCoA2,color=group2,shape=group))+
  geom_point()+theme_bw()+coord_fixed()+
  scale_color_manual(values = c('#AD002A','#42B540'))+
  xlab(temp$anno[1])+ylab(temp$anno[2])
dev.off()

result <- data.frame(fea= c('PCoA1_VE','PCoA2_VE','PCoA1_Juice','PCoA2_Juice'),P1=NA,P2=NA,P3=NA)

k=1
for(i in c('VE','Juice')){
  Info$value <- df$PCoA1
  subinfo <- Info[Info$group2_name==i,]
  result$P1[k] <- round(get_p(subinfo[subinfo$phase2!='Phase2',])$p.value,3)
  result$P2[k] <- round(get_p(subinfo[subinfo$phase2!='Phase0',])$p.value,3)
  result$P3[k] <- round(get_p(subinfo[subinfo$phase2!='Phase1',])$p.value,3)
  k=k+1
  Info$value <- df$PCoA2
  subinfo <- Info[Info$group2_name==i,] #cbw
  result$P1[k] <- round(get_p(subinfo[subinfo$phase2!='Phase2',])$p.value,3)
  result$P2[k] <- round(get_p(subinfo[subinfo$phase2!='Phase0',])$p.value,3)
  result$P3[k] <- round(get_p(subinfo[subinfo$phase2!='Phase1',])$p.value,3)
  k=k+1
}
write.csv(result,'PcoA.csv',quote = F,row.names = F)

k=1
for(i in c('VE','Juice')){
  for(j in paste0('Phase',0:2)){
    d1 <- vegdist(Dat[which(Info$group2_name==i & Info$phase2==j),])
    d1 <- as.numeric(d1)
    assign(paste0('v',k),d1)
    df <- data.frame(group=i,phase=j,value=d1)
    if(k==1){result <- df}
    if(k>1){result <- rbind(result,df)}
    k=k+1
  }
}
wilcox.test(v1,v2);wilcox.test(v1,v3);wilcox.test(v1,v4);wilcox.test(v1,v5);wilcox.test(v1,v6)
wilcox.test(v2,v3);wilcox.test(v2,v4);wilcox.test(v2,v5);wilcox.test(v2,v6)
wilcox.test(v3,v4);wilcox.test(v3,v5);wilcox.test(v3,v6)
wilcox.test(v4,v5);wilcox.test(v4,v6);wilcox.test(v5,v6)

result$comb <-paste0(result$group,result$phase) 

pdf('VJ_bary_boxplot.pdf',width = 6,height = 8)
ggplot(result,aes(x=comb,y=value,fill=phase))+geom_boxplot(outlier.colour = 'grey')+
  theme_bw()+ylab('Bray-Curtis distance')+xlab('')+
  scale_fill_manual(values= c('#4DBBD5','#8491B4','#3C5488'))+
  geom_vline(aes(xintercept = 3.5),linetype = 'dashed',color='grey')
dev.off()

#--------------------------permanova---------------
site <- which(Info$phase2=='Phase0')
adonis2(Dat[site,]~age,data = Info[site,])
adonis2(Dat[site,]~sex,data = Info[site,])

#-----------------------shannon--------------------
kill_0 <- function(X,way=2){
#way = 1 row; way = 2 col
site <- which(apply(X,way,function(x){all(x==0)}))
if(way==1){
if(length(site)>0){X <- X[-site,]}
return(X)
}
if(length(site)>0){X <- X[,-site]}
return(X)
}
Dat <- kill_0(Dat)
sum(row.names(Dat)==Info$seq2)
Info$species_count[1]
sum(Dat[1,]>0)

Info$strain_shannon <- diversity(Dat) 
Info$strain_count <- specnumber(Dat)

#---------------------------------Wilcox test-----------------------------------

Dat <- kill_0(Dat)
fea <- colnames(Dat)
result <- data.frame(group=NA,fea=c(fea,fea),P1=NA,P2=NA,P3=NA,F1=NA,F2=NA,F3=NA)
k=1
for(i in c('VE','Juice')){
  for(j in fea){
    Info$value <- Dat[,j]
    subinfo <- Info[Info$group2_name==i,]
    result$group[k] <- i #cwb
    result$P1[k] <- round(get_p(subinfo[subinfo$phase2!='Phase2',])$p.value,3)
    result$P2[k] <- round(get_p(subinfo[subinfo$phase2!='Phase0',])$p.value,3)
    result$P3[k] <- round(get_p(subinfo[subinfo$phase2!='Phase1',])$p.value,3)
    result$F1[k] <- round(get_FC(subinfo[subinfo$phase2!='Phase2',]),3)
    result$F2[k] <- round(get_FC(subinfo[subinfo$phase2!='Phase0',]),3)
    result$F3[k] <- round(get_FC(subinfo[subinfo$phase2!='Phase1',]),3)
    k=k+1
  }
}

exdata <- read.delim("UHGG_v2_anno.tsv")
exdata <- exdata[exdata$Genome_accession %in% colnames(Dat),]
my_str1 <- function(x,a='g__',b=';s__'){
  y <- unlist(strsplit(x,a))
  y <- y[2]
  z <- unlist(strsplit(y,b))
  return(z[1])
}
exdata$genus <- apply(data.frame(exdata$Lineage),1,my_str1)
my_str1 <- function(x,a='s__',b=';s__'){
  y <- unlist(strsplit(x,a))
  y <- y[2]
  return(y)
}
exdata$species <- apply(data.frame(exdata$Lineage),1,my_str1)
saveRDS(exdata,'micro_anno_v2.rds')

microbe_anno <- readRDS('micro_anno_v2.rds')
df <- result
df[is.na(df)] <- 10
df <- df[df$P1<0.1,]
df$anno <- microbe_anno$Lineage[fsite(df$fea,microbe_anno,'Genome_accession')]
df$species <- microbe_anno$species[fsite(df$fea,microbe_anno,'Genome_accession')]
df$genus <- microbe_anno$genus[fsite(df$fea,microbe_anno,'Genome_accession')]
table(df$group)
table(df$genus[df$group=='Juice'],useNA = 'if')
table(df$genus[df$group=='VE'],useNA = 'if')
write.csv(df,'VJ_meta_wilcox_st.csv',quote = F,row.names = F)
library(ggplot2)

fea <- unique(df$fea)
fea_names <- microbe_anno$label3[fsite(fea,microbe_anno,1)]
n1=7;n2=7
pdf('boxplot significant strains.pdf',width = 4*n1,height = 3*n2)
grid::grid.newpage() 
grid::pushViewport(grid::viewport(layout = grid::grid.layout(n2,n1)))
s1 <- s2 <- 1
for(i in 1:length(fea)){
  Info$value <- Dat[,fea[i]]
  p <- ggplot(Info)+
    geom_boxplot(aes(phase,log(value,10),fill=phase),outlier.color = 'white',show.legend =F,size=0.2)+
    scale_fill_manual(values= c('#4DBBD5','#8491B4','#3C5488'))+
    facet_wrap(~group2_name,ncol =2)+
    theme_classic()+ylab(fea_names[i])+xlab('')+
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


#------------------------------pheatmap-------------------------
library(pheatmap)
exdata <- read.csv("VJ_meta_wilcox_st.csv")
subdata <- Dat[,unique(exdata$fea)]
newdata <- readRDS('ME47_40.rds');temp <- readRDS('MEINFO47.rds')
sum(row.names(newdata)==temp$ID)
temp$value <- newdata$Health.低密度脂蛋白.mmol.L.
newdata <- temp
newdata$id3 <- as.numeric(newdata$id2)
fea <- intersect(newdata$id3,Info$id2)

site <- fsite(fea,Info,'id2')
subdata <- subdata[site,]
row.names(subdata) <- paste0('H',fea)

df <- newdata[,'value',drop=F]
library(psych)
site <- which(newdata$group1=='Juice')
wl <- corr.test(subdata[site,],df[site,],method = 'spearman')
result <- data.frame(J_rho=wl$r,J_p=wl$p)
wl <- corr.test(subdata[-site,],df[-site,],method = 'spearman')
temp <- data.frame(V_rho=wl$r,V_p=wl$p)
result <- cbind(result,temp)
microbe_anno <- readRDS('micro_anno_v2.rds')
fea <- row.names(result)
result$anno <- microbe_anno$Lineage[fsite(fea,microbe_anno,'Genome_accession')]
result$species <- microbe_anno$species[fsite(fea,microbe_anno,'Genome_accession')]
result$genus <- microbe_anno$genus[fsite(fea,microbe_anno,'Genome_accession')]

#---------------------------------boxplot---------------------------------
df <- read.csv("VJ_meta_wilcox_st.csv")
df$species[df$fea=='GUT_GENOME014039'] <- 'Faecalibacterium sp014039'
df$species[df$fea=='GUT_GENOME154740'] <- 'Faecalibacterium sp154740'
subdf <- unique(df[,c(2,10)])
fea <- df$fea[df$group=='Juice']
subdata <- Dat[,fea]
colnames(subdata) <- subdf$species[fsite(fea,subdf,'fea')]
subdata$id <- row.names(subdata)
result <- melt(subdata)
colnames(result) <- c('id','fea','value')
result$phase <- Info$phase2[fsite(result$id,Info,'seq2')]
result$value2 <- log10(result$value)
sub2data <- subdata[-ncol(subdata)]
value <- colSums(sub2data[Info$phase2=='Phase1',])
fea <- colnames(sub2data)[order(value,decreasing = T)][1:15]
fea
subresult <- result[fsite(fea,result,'fea'),]
sub2data <- sub2data[,fea]
v1 <- colMeans(sub2data)
subresult$fea <- factor(subresult$fea,levels = fea[order(v1,decreasing = T)])

pdf('Juice top 15 strain.pdf',width = 7,height = 4.5)
ggplot(subresult, aes(x=fea, y=value2, fill=phase)) + 
  geom_boxplot(outlier.colour = 'grey',outlier.size = 0.2)+
  scale_fill_manual(values= c('#4DBBD5','#8491B4','#3C5488'))+
  theme_classic()+xlab('')+ylab('Relative abundance (log10)')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

fea <- df$fea[df$group=='VE']
subdata <- Dat[,fea]
colnames(subdata) <- subdf$species[fsite(fea,subdf,'fea')]
subdata$id <- row.names(subdata)
result <- melt(subdata)
colnames(result) <- c('id','fea','value')
result$phase <- Info$phase2[fsite(result$id,Info,'seq2')]
result$value2 <- log10(result$value)
v1 <- colMeans(subdata[-ncol(subdata)])
fea <- colnames(subdata)[1:10]
result$fea <- factor(result$fea,levels = fea[order(v1,decreasing = T)])
pdf('VE top 10 strain.pdf',width = 6,height = 4.5)
ggplot(result, aes(x=fea, y=value2, fill=phase)) + 
  geom_boxplot(outlier.colour = 'grey',outlier.size = 0.2)+
  scale_fill_manual(values= c('#4DBBD5','#8491B4','#3C5488'))+
  theme_classic()+xlab('')+ylab('Relative abundance (log10)')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()