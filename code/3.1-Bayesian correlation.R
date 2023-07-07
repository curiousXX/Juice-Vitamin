source('F_know_data.R')
setwd('YOUR_PATH')
topx <- function(value,x=3){
  y <- order(value)
  s1 <- y[1:3]
  y <- order(value,decreasing = T)
  s2 <- y[1:3]
  return(c(s1,s2))
}
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

Info <- readRDS('MEINFO47.rds')
Dat <- readRDS('ME47_45.rds')
microbe_anno <- readRDS('micro_anno_v2.rds')
Info$value <- Dat$Health.低密度脂蛋白.mmol.L.


library(BayesFactor)
v1 <- value
v2 <- Info$value
ggplot()+geom_line(mapping = aes(x=v2,y=v1),color='black')+
  geom_line(mapping = aes(x=v2,y=log10(v1+1)),color='red')+
  geom_line(mapping = aes(x=v2,y=log10(v1+1e-10)),color='blue')+theme_classic()
#-------------------------VE---------------------------
exdata <- read.csv("VE_LDL_change_ratio.csv")
temp <- Info[Info$phase2!=2,]
id <- temp$seq2[fsite(exdata$names,temp,'name3')]
Dat <- readRDS('D112_1098.rds')
Dat <- Dat[id,]
Dat <- kill_0(Dat)
subinfo <- Info[fsite(row.names(Dat),Info,'seq2'),]

site <- c()
s1 <- which(subinfo$phase2==0)
for(i in 1:ncol(Dat)){
  value <- Dat[,i]
  value <- as.numeric(value>0)
  if(sum(value[s1])<5 | sum(value[-s1])<5){site <- c(site,i)}
}
Dat <- Dat[,-site]
saveRDS(Dat,'VE_jq.rds')
exdata <- read.csv("VJ_meta_wilcox_st.csv")
Dat <- Dat[,unique(exdata$fea[exdata$group=='VE'])]
result <- data.frame(features=colnames(Dat),rho1=NA,rho1_95 = NA, p1_strong=NA,p1_overzero=NA,
                     rho2=NA,rho2_95 = NA,p2_strong=NA,p2_overzero=NA)
fea <- result$features
saveRDS(subinfo,'VE_info_jq.rds')
for(i in 1:length(fea)){
  print(paste0('-------------------i = ',i,'---------------------'))
  subinfo$value2 <- Dat[,fea[i]] 
  v1 <- subinfo$value[subinfo$phase2==0]
  v2 <- subinfo$value[subinfo$phase2==1]
  v3 <- subinfo$value2[subinfo$phase2==0]
  v4 <- subinfo$value2[subinfo$phase2==1]
  v3 <- log10(v3+1e-10) #cbw
  v4 <- log10(v4+1e-10)
 
  fit = correlationBF(x=v1,y=v3,posterior = TRUE, iterations = 5000)
  value <- fit[,"rho"]
  hpd95 = HPDinterval(as.mcmc(as.numeric(value)), prob=0.95)

  result$rho1[i] <- mean(value)
  result$rho1_95[i] <- paste0(round(hpd95[,"lower"],3),'~',round(hpd95[,"upper"],3))
  result$p1_strong[i] <- 1-mean(value > -0.1 & value < 0.1)
  result$p1_overzero[i] <- mean(value >= 0)
  
  fit = correlationBF(x=v2,y=v4,posterior = TRUE, iterations = 5000)
  value <- fit[,"rho"]
  hpd95 = HPDinterval(as.mcmc(as.numeric(value)), prob=0.95)
  
  result$rho2[i] <- mean(value)
  result$rho2_95[i] <- paste0(round(hpd95[,"lower"],3),'~',round(hpd95[,"upper"],3))
  result$p2_strong[i] <- 1-mean(value > -0.1 & value < 0.1)
  result$p2_overzero[i] <- mean(value >= 0)
}

result <- data.frame(features=colnames(Dat),rho1=NA, p1_strong=NA,p1_overzero=NA,
                     rho2=NA,p2_strong=NA,p2_overzero=NA)
for(i in 1:length(fea)){
  print(paste0('-------------------i = ',i,'---------------------'))
  subinfo$value2 <- Dat[,fea[i]] 
  v1 <- subinfo$value[subinfo$phase2==0]
  v2 <- subinfo$value[subinfo$phase2==1]
  v3 <- subinfo$value2[subinfo$phase2==0]
  v4 <- subinfo$value2[subinfo$phase2==1]
  
  wp <- rob.cor.mcmc(x=data.frame(LDL=v1,Strain=v3),iter = 2000,warmup = 1000,chains=4,treedepth = 10)
  result$rho1[i] <- wp$rho
  result$p1_strong[i] <- 1-wp$p_weak
  result$p1_overzero[i] <- wp$p_overzero 
  
  wp <- rob.cor.mcmc(x=data.frame(LDL=v2,Strain=v4),iter = 2000,warmup = 1000,chains=4,treedepth = 10)
  result$rho2[i] <- wp$rho
  result$p2_strong[i] <- 1-wp$p_weak
  result$p2_overzero[i] <- wp$p_overzero
  
}
result$species<- microbe_anno$species[fsite(result$features,microbe_anno,13)]

write.csv(result,'VE_LDL_strain_rho_BC_signi.csv',quote = F,row.names = F)

result <- read.csv("VE_LDL_strain_rho_BC.csv")
subresult <- result
subresult[is.na(subresult)] <- 10
site <- which(subresult$p1_overzero<0.1 | subresult$p2_overzero<0.1 |
                subresult$p1_overzero>0.9 | subresult$p2_overzero>0.9)
subresult <- subresult[site,]
fea <- subresult$features
site <- which(subinfo$phase2==0)

n1=2;n2=length(fea)
pdf('VE BC points plot.pdf',width = 3*n1,height = 3*n2)
grid::grid.newpage() 
grid::pushViewport(grid::viewport(layout = grid::grid.layout(n2,n1)))
s1 <- s2 <- 1
for(i in 1:length(fea)){
  subinfo$value2 <- Dat[,fea[i]]
  df = subinfo[site,] #cbw
  v1 <- round(subresult$rho1[i],3)
  v2 <- round(subresult$p1_overzero[i],3)
  p <- ggplot(data = df,mapping = aes(x=value,y=log10(value2+1e-10)))+
                geom_point()+theme_bw()+xlab('LDL-C (mmol/L)')+
                ylab(paste0(subresult$species[i],' (log10)'))+
                labs(title = paste0('r = ',v1,' p = ',v2))
  print(p,vp = vplayout(s1,s2))
  if(s2==n1){
    s1 <- s1+1
    s2=1
  }
  else{s2=s2+1} 
  
  df = subinfo[-site,]
  v1 <- round(subresult$rho2[i],3)
  v2 <- round(subresult$p2_overzero[i],3)
  p <- ggplot(data = df,mapping = aes(x=value,y=log10(value2+1e-10)))+
    geom_point()+theme_bw()+xlab('LDL-C (mmol/L)')+
    ylab(paste0(subresult$species[i],' (log10)'))+
    labs(title = paste0('r = ',v1,' p = ',v2))
  print(p,vp = vplayout(s1,s2))
  if(s2==n1){
    s1 <- s1+1
    s2=1
  }
  else{s2=s2+1} 
  
}
dev.off()

subdata <- Dat[,subresult$features]
library(pheatmap)
min(subdata[subdata>0])
site <- order(subinfo$phase2)
subdata <- subdata[site,];subinfo <- subinfo[site,]
sum(row.names(subdata)==subinfo$seq2)
colnames(subdata) <- microbe_anno$species[fsite(colnames(subdata),microbe_anno,13)]

anno_row <- data.frame(stage=factor(subinfo$phase2),row.names = row.names(subdata))
pheatmap(log10(subdata+1e-10),show_rownames = F,
         annotation_row = anno_row,cluster_rows = F)

subresult$species <- microbe_anno$species[fsite(subresult$features,microbe_anno,13)]
df <- subresult[,c(2,4)];row.names(df) <- subresult$species
df2 <- subresult[,c(3,5)]
temp <- df2
df2[temp<0.1] <- '.'
df2[temp<0.05] <- '*'
df2[temp<0.01] <- '**'
df2[temp>=0.1] <- ''
colnames(df) <- c('Baseline','Intervention')
pheatmap(df,display_numbers = df2,cluster_cols = F,
         angle_col = 0,cellwidth = 39,height = 10,
         fontsize = 8)
exdata <- read.csv("VJ_meta_wilcox_st.csv")
fea <- exdata$species[exdata$group=='VE']
site <- which(row.names(df) %in% fea)
intersect(row.names(df),fea)


#------------------------Juice-------------------------
exdata <- read.csv("Juice_LDL_change_ratio.csv")
temp <- Info[Info$phase2!=2,]
id <- temp$seq2[fsite(exdata$names,temp,'name3')]
Dat <- readRDS('D112_1098.rds')
Dat <- Dat[id,]
Dat <- kill_0(Dat)
site <- c()

for(i in 1:ncol(Dat)){
  value <- Dat[,i]
  value <- as.numeric(value>0)
  if(sum(value)<5){site <- c(site,i)}
}
Dat <- Dat[,-site]

fea <- result$features
subinfo <- Info[fsite(row.names(Dat),Info,'seq2'),]

saveRDS(Dat,'Juice_jq.rds')
saveRDS(subinfo,'Juice_info_jq.rds')

result <- data.frame(features=colnames(Dat),rho1=NA,rho1_95 = NA, p1_strong=NA,p1_overzero=NA,
                     rho2=NA,rho2_95 = NA,p2_strong=NA,p2_overzero=NA,
                     rho3=NA,rho3_95 = NA,p3_strong=NA,p3_overzero=NA)
fea <- result$features
for(i in 1:length(fea)){
  print(paste0('-------------------i = ',i,'---------------------'))
  subinfo$value2 <- Dat[,fea[i]] 
  v1 <- subinfo$value[subinfo$phase2==0]
  v2 <- subinfo$value[subinfo$phase2==1]
  v3 <- subinfo$value2[subinfo$phase2==0]
  v4 <- subinfo$value2[subinfo$phase2==1] 
  v5 <- v2-v1
  v6 <- v4-v3
  v3 <- log10(v3+1e-10) #cwb
  v4 <- log10(v4+1e-10)

  fit = correlationBF(x=v1,y=v3,posterior = TRUE, iterations = 5000)
  value <- fit[,"rho"]
  hpd95 = HPDinterval(as.mcmc(as.numeric(value)), prob=0.95)
  
  result$rho1[i] <- mean(value)
  result$rho1_95[i] <- paste0(round(hpd95[,"lower"],3),'~',round(hpd95[,"upper"],3))
  result$p1_strong[i] <- 1-mean(value > -0.1 & value < 0.1)
  result$p1_overzero[i] <- mean(value >= 0)
  
  fit = correlationBF(x=v2,y=v4,posterior = TRUE, iterations = 5000)
  value <- fit[,"rho"]
  hpd95 = HPDinterval(as.mcmc(as.numeric(value)), prob=0.95)
  
  result$rho2[i] <- mean(value)
  result$rho2_95[i] <- paste0(round(hpd95[,"lower"],3),'~',round(hpd95[,"upper"],3))
  result$p2_strong[i] <- 1-mean(value > -0.1 & value < 0.1)
  result$p2_overzero[i] <- mean(value >= 0)
  
  fit = correlationBF(x=v5,y=v6,posterior = TRUE, iterations = 5000)
  value <- fit[,"rho"]
  hpd95 = HPDinterval(as.mcmc(as.numeric(value)), prob=0.95)
  
  result$rho3[i] <- mean(value)
  result$rho3_95[i] <- paste0(round(hpd95[,"lower"],3),'~',round(hpd95[,"upper"],3))
  result$p3_strong[i] <- 1-mean(value > -0.1 & value < 0.1)
  result$p3_overzero[i] <- mean(value >= 0)
  
}

result$species <- microbe_anno$species[fsite(result$features,microbe_anno,13)]
write.csv(result,'Juice_LDL_strain_BC.csv',quote = F,row.names = F)

site <- c()
s1 <- which(subinfo$phase2==0)
for(i in 1:ncol(Dat)){
  value <- Dat[,i]
  value <- as.numeric(value>0)
  v1 <- value[s1]
  v2 <- value[-s1]
  if(sum(v1)<5 & sum(v2)<5){site <- c(site,i)}
}
fea <- colnames(Dat)[site]
subresult <- result[-fsite(fea,result,'features'),]
site <- which(subresult$p1_overzero<0.1 | subresult$p2_overzero<0.1 |
                subresult$p1_overzero>0.9 | subresult$p2_overzero>0.9)
subresult <- subresult[site,]
site <- which(subresult$rho1 * subresult$rho2 >0)
subresult <- subresult[site,]
subresult <- subresult[which(abs(subresult$rho1)>0.1 & abs(subresult$rho2)>0.1),]
fea <- subresult$features
site <- which(subinfo$phase2==0)

n1=4;n2=ceiling(length(fea)/2)
pdf('Juice BC points plot.pdf',width = 3*n1,height = 3*n2)
grid::grid.newpage() 
grid::pushViewport(grid::viewport(layout = grid::grid.layout(n2,n1)))
s1 <- s2 <- 1
for(i in 1:length(fea)){
  subinfo$value2 <- Dat[,fea[i]]
  df = subinfo[site,]
  v1 <- round(subresult$rho1[i],3)
  v2 <- round(subresult$p1_overzero[i],3)
  p <- ggplot(data = df,mapping = aes(x=value,y=log10(value2+1e-10)))+
    geom_point()+theme_bw()+xlab('LDL-C (mmol/L)')+
    ylab(paste0(subresult$species[i],' (log10)'))+
    labs(title = paste0('r = ',v1,' p = ',v2))
  print(p,vp = vplayout(s1,s2))
  if(s2==n1){
    s1 <- s1+1
    s2=1
  }
  else{s2=s2+1} 
  
  df = subinfo[-site,]
  v1 <- round(subresult$rho2[i],3)
  v2 <- round(subresult$p2_overzero[i],3)
  p <- ggplot(data = df,mapping = aes(x=value,y=log10(value2+1e-10)))+
    geom_point()+theme_bw()+xlab('LDL-C (mmol/L)')+
    ylab(paste0(subresult$species[i],' (log10)'))+
    labs(title = paste0('r = ',v1,' p = ',v2))
  print(p,vp = vplayout(s1,s2))
  if(s2==n1){
    s1 <- s1+1
    s2=1
  }
  else{s2=s2+1} 
  
}
dev.off()

subresult <- result
site <- which(subresult$p3_overzero<0.1  | subresult$p3_overzero>0.9 )
subresult <- subresult[site,]
fea <- subresult$features
n1=4;n2=ceiling(length(fea)/2)
pdf('Juice change BC points plot.pdf',width = 3*n1,height = 3*n2)
grid::grid.newpage() 
grid::pushViewport(grid::viewport(layout = grid::grid.layout(n2,n1)))
s1 <- s2 <- 1
for(i in 1:length(fea)){
  subinfo$value2 <- Dat[,fea[i]]
  
  df = subinfo[subinfo$value2>0,]
  df$phase2 <- as.factor(df$phase2)
  v1 <- round(subresult$rho3[i],3)
  v2 <- round(subresult$p3_overzero[i],3)
  p <- ggplot(data = df,mapping = aes(x=value,y=log10(value2),color=phase2))+
    geom_point()+theme_bw()+xlab('LDL-C (mmol/L)')+
    ylab(paste0(subresult$species[i],' (log10)'))+
    labs(title = paste0('r = ',v1,' p = ',v2))
  print(p,vp = vplayout(s1,s2))
  if(s2==n1){
    s1 <- s1+1
    s2=1
  }
  else{s2=s2+1} 
  
}
dev.off()


subdata <- Dat[,subresult$features]
exdata <- read.csv("VJ_meta_wilcox_st.csv")
fea <- exdata$species[exdata$group=='Juice']
intersect(subresult$features,fea)



