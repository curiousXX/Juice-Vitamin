#---------------------------------说明--------------------------------------
#该函数可以对输入数据的特征（列变量）做划分，得到连续变量，二分类或单分类变量、多分类变量
#factor_action主要用于区分是连续变量还是分类变量，如果您的数据中有某个分类变量是20分类，请设置factor_action在20或以上
#如果希望检查数据的NA情况，请设置 NA_action='extra'，将输出所有特征的NA数量
#如果希望更直观地查看数据，可以通过 draw=TRUE进行画图，图片将直接以pdf输出，如果您希望确定文件名，请设置参数File
#本函数运行6200*112数据仅需要21秒,如果运行时间过长请检查输入数据是否出错

#---------------------------主函数------------------------------------------
know_data <- function(X,factor_action=10,NA_action='ignore',File='auto',draw=FALSE,silent=FALSE){
  time1<-Sys.time()
  #检查R包
  require(ggplot2)
  require(grid)
  require(nortest)
  require(dplyr)
  #检查参数
  if(class(X)=='matrix'){X <- as.data.frame(X)}
  if(class(X)!='data.frame'){return('Please input matrix or data frame.')}
  if(!silent){print('factor_action : 10(default) or a number.')}
  if(class(factor_action)!='numeric'){return('Please input correct parm(factor_action).')}
  
  if(!silent){print('NA_action : ignore(default), extra')}
  if(NA_action!='ignore' & NA_action!='extra'){return('Please input correct parm(NA_action).')}
  
  if(!silent){print('File : auto(default) or filename.')}
  if(File!='auto' & class(File)!='character'){return('Please input correct parm(File).')}
  
  if(!silent){print('draw : FALSE(default) or TRUE')}
  if(!is.logical(draw)){return('Please input correct parm(draw).')}
  #第一步 ： 确定变量性质，离散 or 连续 ，得到列表
  fea_result <- feature_split(X,pd=factor_action)
  #第二步：画图，连续变量画密度图，离散变量画饼图,输出pdf
  if(draw){
    Pname1 <- paste0(File,'_Contin.pdf')
    Pname2 <- paste0(File,'_Discrete.pdf')
    if(File=='auto'){
      Pname1 <- paste0('konw_data','_Contin.pdf')
      Pname2 <- paste0('konw_data','_Discrete.pdf')
    }
    Cfea <- as.character(fea_result$Contin_feature)
    Dfea <- c(fea_result$`Discrete_feature(level<=2)`,fea_result$`Discrete_feature(level>2)`)
    Dfea <- as.character(Dfea)
    ifelse(nrow(X)>2000, Method <- 'lillie' , Method <- 'shapiro')
    site <- fsite(Cfea,X,name = 'col')
    if(length(site)>0){
      h <- w <- 5*(ceiling(sqrt(length(site))))
      pdf(file = Pname1,height = h,width = w,family = 'GB1')
      Gaussian_plot(X,site = site,Method = Method)
      dev.off()
    }
    site <- fsite(Dfea,X,name = 'col')
    if(length(site)>0){
      h <- w <- 10*(ceiling(sqrt(length(site))))
      pdf(file = Pname2,height = h,width = w,family = 'Helvetica')
      Agur_plot(X,site=site)
      dev.off()
    }
  }
  #第三步：输出
  if(NA_action=='extra'){
    na_result <- na_count(X)
    result <- list(fea_result,na_result)
    names(result) <- c('fea_result','na_result')
    time2<-Sys.time()
    if(!silent){print(paste0('Runing time : ',(time2-time1)))}
    return(result)
  }
  time2<-Sys.time()
  if(!silent){print(paste0('Runing time : ',(time2-time1)))}
  return(fea_result)
}

#--------------------子函数--------------------------------------------
#辅助小函数
#众数
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#确定位置
fsite <- function(v,m,colnum=NULL,name="row"){
  v <- as.character(v)
  if(is.character(colnum)){colnum <- which(colnames(m)==colnum)}
  if(is.numeric(colnum)){
    site.name=m[,colnum]
    site=c()
    for(i in 1:length(v)){
      x=which(site.name==v[i])
      site=c(site,x)
    }
    return(site)
  }
  if(name!="row" && name!= "col"){
    print("error:parameter name need to be row or col")
  }
  else{
    if(name=="row"){site.name=row.names(m)}
    if(name=="col"){site.name=colnames(m)}
    site=c()
    for(i in 1:length(v)){
      x=which(site.name==v[i])
      site=c(site,x)
    }
    return(site)
  }
}

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols),byrow =TRUE)
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp =grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


vplayout <- function(x,y){grid::viewport(layout.pos.row = x, layout.pos.col = y)}  
na_count <- function(Dat,result=TRUE,view = FALSE){
  print(paste0('NA ratio : ',sum(is.na(Dat))/(ncol(Dat)*nrow(Dat))))
  #-----row na count-----------
  df_row <- apply(Dat,1,function(x){sum(is.na(x))})
  df_row <- as.data.frame(df_row)
  colnames(df_row) <- 'count'
  #df_row$name=row.names(df_row)
  df_row$count <- as.numeric(as.character(df_row$count))
  df_row$perc <- round((df_row$count)/ncol(Dat),2)
  df_row <- df_row[order(df_row$count,decreasing = TRUE),]
  
  #-----col na count-----------
  df_col <- apply(Dat,2,function(x){sum(is.na(x))})
  df_col <- as.data.frame(df_col)
  colnames(df_col) <- 'count'
  #df_col$name=row.names(df_col)
  df_col$count <- as.numeric(as.character(df_col$count))
  df_col$perc <- round((df_col$count)/nrow(Dat),2)
  df_col <- df_col[order(df_col$count,decreasing = TRUE),]
  
  #-------
  if(view){
    View(df_row)
    View(df_col)
  }
  if(result){
    result_output <- list(df_row,df_col)
    names(result_output) <- c('row','col')
    return(result_output)
  }
  return(NULL)
}
#---------------------------第一步----------------------------------
feature_split <- function(Dat,pd=5){
  C.check <- function(x){
    x <- as.numeric(na.omit(as.character(x)))
    count=sum(is.na(x))
    if(count<(length(x)*0.5) & length(unique(x))>pd){return(1)}
    return(0)
  }
  is_contin <- as.data.frame(apply(Dat,2,C.check))
  colnames(is_contin) <- 'result'
  site=which(is_contin$result==1)
  ifelse(length(site)>0, contin_fea <- colnames(Dat)[site], contin_fea <- '')
  
  
  site=which(is_contin$result==0)
  if(length(site)>0){discrete_Dat=Dat[,site]}
  if(length(site)==0){
    discrete_fea_two <- discrete_fea_more <- ''
    result <- list(contin_fea,discrete_fea_two,discrete_fea_more)
    names(result) <- c('Contin_feature','Discrete_feature(level<=2)','Discrete_feature(level>2)')
    return(result)
  }
  str2 <- function(x){
    x <- na.omit(as.character(x))
    if(length(levels(factor(x)))<=2){return(1)}
    return(0)}
    
  is_two <- as.data.frame(apply(discrete_Dat,2,str2))
  colnames(is_two) <- 'result'
  site=fsite(row.names(is_two)[is_two$result==1],discrete_Dat,name='col')
  ifelse(length(site)>0, discrete_fea_two<- colnames(discrete_Dat)[site], discrete_fea_two <- '')
    
  site=fsite(row.names(is_two)[is_two$result==0],discrete_Dat,name='col')
  ifelse(length(site)>0, discrete_fea_more<- colnames(discrete_Dat)[site], discrete_fea_more <- '')

  result <- list(contin_fea,discrete_fea_two,discrete_fea_more)
  names(result) <- c('Contin_feature','Discrete_feature(level<=2)','Discrete_feature(level>2)')
  return(result)
}
#---------------------------------第二步-------------------------------------
Gaussian_plot  <- function(data,site=NULL,Method='lillie',beauty=FALSE){
  cnames <- colnames(data)
  if(!is.null(site)){
    cnames=cnames[site]
  }
  n=ceiling(sqrt(length(cnames)))
  #if(n*n==length(cnames)){m=n}
  grid::grid.newpage()  ###新建图表版面
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n,n))) #将版面分成2*2矩阵
  
  s1 <- s2 <- 1
  #离群值判断：大于中位数加上四分位的值
  for(i in 1:length(cnames)){
    Check_point <- 0
    value <- as.numeric(as.character(data[,cnames[i]]))
    data[,cnames[i]] <- value
    temp <- na.omit(value)
    if(max(temp)>(10*median(temp)) & median(temp)>0){
      print(paste0(cnames[i],' has outlier'))
      xlim_r=median(temp)+quantile(temp)[4]
      xlim_l=min(temp)
      Check_point <- 1
    }
    
    if(Method=='shapiro'){p_value <- shapiro.test(value)$p.value}
    if(Method=='lillie'){p_value <- (nortest::lillie.test(value))[["p.value"]]}
    p_value <- signif(p_value)
    p=ggplot2::ggplot(data,ggplot2::aes_string(x=cnames[i]))+ggplot2::geom_line(stat = "density")+  
      ggplot2::labs(title = paste0('p=',p_value))+theme_classic()
    if(Check_point ==1 & beauty==TRUE){p=p+ggplot2::xlim(xlim_l,xlim_r)}
    print(p,vp = vplayout(s1,s2))
    if(s2==n){
      s1 <- s1+1
      s2=1
    }
    else{s2=s2+1}
  }
  #return(p_result)
}

Agur_plot <- function(Data,site=NULL){
  cnames <- colnames(Data)
  if(!is.null(site)){
    cnames=cnames[site]
  }
  n=ceiling(sqrt(length(cnames)))
  n=ceiling(sqrt(length(cnames)))
  grid::grid.newpage()  ###新建图表版面
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n,n))) #将版面分成2*2矩阵
  
  s1=s2=1
  for(i in 1:length(cnames)){
    
    value <- Data[,cnames[i]]
    tvalue <-table(value) 
    df <- data.frame(value = as.numeric(tvalue),
                     Group = names(tvalue) )%>%
      # factor levels need to be the opposite order of the cumulative sum of the values
      mutate(Group = factor(Group),
             cumulative = cumsum(value),
             midpoint = cumulative - value / 2,
             label = paste0(Group, " ", round(value / sum(value) * 100, 2), "% (",value,')'))
    
    p <- ggplot(df, aes(x = 1, weight = value, fill = Group)) +
      geom_bar(width = 1, position = "stack") +
      coord_polar(theta = "y")  +
      labs(x = '', y = cnames[i], title = '') + 
      theme(legend.title = element_blank(), legend.position = "bottom")+
      scale_fill_discrete(breaks = df$Group, labels = df$label)+
      theme(axis.ticks = element_blank(),axis.text.x = element_blank(),
            axis.text.y = element_blank()) + 
      theme(panel.grid=element_blank(),panel.border=element_blank())
    print(p,vp = vplayout(s1,s2))
    if(s2==n){
      s1 <- s1+1
      s2=1
    }
    else{s2=s2+1}
  }
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
set.seed(21)
library(ggplot2)
print('Already set seed: 21')
print('Already load ggplot2')
print('Welcome back, Victor.')
