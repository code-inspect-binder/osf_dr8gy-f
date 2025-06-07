SG<-function(x){
  if(is.integer(DATA[,match(x,names(DATA))])){
  
  a<-quantile(DATA[,match(x,names(DATA))],na.rm=T);
  
  b<-paste(x,sep = "",">",a[3]) ;
  
  c<-DATA[DATA[x]>a[3],];
  
  d <- twoby2( 1- as.numeric(as.character(c$CGR)), 1- as.numeric(as.character(c$status_1an)));
  
  e<-ifelse(d$p.value[2]<0.0001, "p<0.0001", round(d$p.value[2],4));
  
  
  b1<-paste(x,sep = "","<=",a[3]) ;
  
  c1<-DATA[DATA[x]<=a[3],];
  
  d1 <- twoby2( 1- as.numeric(as.character(c1$CGR)), 1- as.numeric(as.character(c1$status_1an)));
  
  e1 <-ifelse(d1$p.value[2]<0.0001, "p<0.0001", round(d1$p.value[2],4));
  
  res<-c(paste(b),
         paste(round(d$measures[2,1],2)),
         paste(round(d$measures[2,2],2)),
         paste(round(d$measures[2,3],2)),
         paste(e));
  res1<-c(paste(b1),
          paste(round(d1$measures[2,1],2)),
          paste(round(d1$measures[2,2],2)),
          paste(round(d1$measures[2,3],2)),
          paste(e1))
  MAT<-matrix(NA,2,5);
  MAT[1,]<-res;
  MAT[2,]<-res1;
  } 
else{
    b2<-paste(x);
  c2<-DATA[DATA[x]==1,];
  d2<- twoby2( 1- as.numeric(as.character(c2$CGR)), 1- as.numeric(as.character(c2$status_1an)));
  
  e2<-ifelse(d2$p.value[2]<0.0001, "p<0.0001", round(d2$p.value[2],4));
  res2<-c(paste(b2),
          paste(round(d2$measures[2,1],2)),
          paste(round(d2$measures[2,2],2)),
          paste(round(d2$measures[2,3],2)),
          paste(e2));
  MAT<-matrix(NA,1,5);
  MAT[1,]<-res2}
  return(MAT);}


SG1<-function(x){
  if(is.integer(PS[,match(x,names(PS))])){
    
    a<-quantile(PS[,match(x,names(PS))],na.rm=T);
    
    b<-paste(x,sep = "",">",a[3]) ;
    
    c<-PS[PS[x]>a[3],];
    
    d <- twoby2( 1- as.numeric(as.character(c$CGR)), 1- as.numeric(as.character(c$status_1an)));
    
    e<-ifelse(d$p.value[2]<0.0001, "p<0.0001", round(d$p.value[2],4));
    
    
    b1<-paste(x,sep = "","<=",a[3]) ;
    
    c1<-PS[PS[x]<=a[3],];
    
    d1 <- twoby2( 1- as.numeric(as.character(c1$CGR)), 1- as.numeric(as.character(c1$status_1an)));
    
    e1 <-ifelse(d1$p.value[2]<0.0001, "p<0.0001", round(d1$p.value[2],4));
    
    res<-c(paste(b),
           paste(round(d$measures[2,1],2)),
           paste(round(d$measures[2,2],2)),
           paste(round(d$measures[2,3],2)),
           paste(e));
    res1<-c(paste(b1),
            paste(round(d1$measures[2,1],2)),
            paste(round(d1$measures[2,2],2)),
            paste(round(d1$measures[2,3],2)),
            paste(e1))
    MAT<-matrix(NA,2,5);
    MAT[1,]<-res;
    MAT[2,]<-res1;
  } 
  else{
    b2<-paste(x);
    c2<-PS[PS[x]==1,];
    d2<- twoby2( 1- as.numeric(as.character(c2$CGR)), 1- as.numeric(as.character(c2$status_1an)));
    
    e2<-ifelse(d2$p.value[2]<0.0001, "p<0.0001", round(d2$p.value[2],4));
    res2<-c(paste(b2),
            paste(round(d2$measures[2,1],2)),
            paste(round(d2$measures[2,2],2)),
            paste(round(d2$measures[2,3],2)),
            paste(e2));
    MAT<-matrix(NA,1,5);
    MAT[1,]<-res2}
  return(MAT);}


GRA<-function(x, y, labY, title, yPospval){
  prop_dat <- na.omit(DATA[,c(x, y)])
  colnames(prop_dat)<-c("x","y")
  WC<-wilcox.test(y~x, data=prop_dat)
  pval<-ifelse(WC$p.value<0.0001,"<0.0001",WC$p.value)
  mods_plot<-ggplot(prop_dat, aes(x=x, y=log(y), fill=x))+
    geom_boxplot() +  
    xlab("")+
    ylab(labY)+
    theme(axis.title.x = element_text( face="bold", size=20),
          axis.title.y = element_text( face="bold", size=20),
          axis.text.x = element_text(face="bold", size=20),
          axis.text.y = element_text(face="bold", size=20 ),
          legend.position = "none")+
    scale_fill_manual(breaks = c("No transfusion", "Transfusion" ), 
                      values=c("blue","red"))+
    annotate("text", x = 1, y = yPospval, label = paste("p=",pval),hjust = 0)
  print(mods_plot)
}

GRAm<-function(x, y, labY, title, yPospval){
  prop_dat <- na.omit(PS[,c(x, y)])
  colnames(prop_dat)<-c("x","y")
  WC<-wilcox.test(y~x, data=prop_dat)
  pval<-ifelse(WC$p.value<0.0001,"<0.0001",round(WC$p.value,4))
  mods_plot<-ggplot(prop_dat, aes(x=x, y=log(y), fill=x))+
    geom_boxplot() +  
    xlab("")+
    ylab(labY)+
    theme(axis.title.x = element_text( face="bold", size=20),
          axis.title.y = element_text( face="bold", size=20),
          axis.text.x = element_text(face="bold", size=20),
          axis.text.y = element_text(face="bold", size=20 ),
          legend.position = "none")+
    scale_fill_manual(breaks = c("Transfusion", "No transfusion"), 
                      values=c("red", "blue"))+
    annotate("text", x = 1, y = yPospval, label = paste("p=",pval),hjust = 0)
  print(mods_plot)
}


INT<-function(x){
  df<-DATAtab1 %>%
    dplyr::select(CGR,status_1an,x);
  require(dplyr);
  if(is.numeric(df[,match(x,names(df))])){
    med<-quantile(df[,match(x,names(df))],na.rm=T)
    df<-df%>%
      mutate(x_cat=
               if_else(df[,match(x,names(df))]<= med[[3]],
                       paste0(x,"<=",med[[3]]),
                       paste0(x,">",med[[3]])));
    df$x_cat<-as.factor(as.character(df$x_cat))                     
  
    mod1<-glm(status_1an ~ CGR *x_cat, data=df,family = binomial());
    #pval interaction
    require(car);
    anov<-Anova(mod1, type=3);
    pval<-round(anov$`Pr(>Chisq)`[3],2);
    
    #OR and estimate first level
    require(gtsummary);
    OR1<-(tbl_regression(mod1, exponentiate = TRUE));
    ci1l<-round(OR1$table_body$conf.low[OR1$table_body$label==paste0("CGR")],2);
    ci1h<-round(OR1$table_body$conf.high[OR1$table_body$label==paste0("CGR")],2);    
    est1<-round(OR1$table_body$estimate[OR1$table_body$label==paste0("CGR")],1);
    
    df$x_cat1<-relevel(df$x_cat, ref= levels(df$x_cat)[2])
    
    mod2<-glm(status_1an ~ CGR * x_cat1, data = df, family = binomial())
    
    OR2<-(tbl_regression(mod2, exponentiate = TRUE))
    ci2l<-round(OR2$table_body$conf.low[OR2$table_body$label==paste0("CGR")],2);
    ci2h<-round(OR2$table_body$conf.high[OR2$table_body$label==paste0("CGR")],2);    
    est2<-round(OR2$table_body$estimate[OR2$table_body$label==paste0("CGR")],2)
    
    
    res1<-c(paste0(levels(df$x_cat)[1]),
           est1,
           ci1l,
           ci1h,
           pval);
    res2<-c(paste0(levels(df$x_cat1)[1]),
            est2,
            ci2l,
            ci2h,
            NA)

    MAT<-matrix(NA,2,5);
    MAT[1,]<-res1;
    MAT[2,]<-res2;
  } 
  
  else{
    df$x_cat2 <-as.factor(as.character(df[,match(x,names(df))])) #[1] "0" "1"
    
    df<-na.omit(df)
    
    mod3<-glm(status_1an ~ CGR * x_cat2, data=df, family = binomial());
    
    #pval interaction
    require(car);
    anov1<-Anova(mod3, type=3);
    pval1<- ifelse(anov1$`Pr(>Chisq)`[3]<0.0001, "p<0.0001", round(anov1$`Pr(>Chisq)`[3],4));
    
    #OR and estimate first level
    require(gtsummary);
    OR3<-(tbl_regression(mod3, exponentiate = TRUE));
    
    ci3l<-round(OR3$table_body$conf.low[OR3$table_body$label==paste0("CGR")],2);
    ci3h<-round(OR3$table_body$conf.high[OR3$table_body$label==paste0("CGR")],2);    
    est3<-round(OR3$table_body$estimate[OR3$table_body$label==paste0("CGR")],2);
    
    #OR and estimate second level
    df$x_cat3<-relevel(df$x_cat2, ref= levels(df$x_cat2)[2]);
    
    mod4<-glm(status_1an ~ CGR * x_cat3, data = df, family = binomial());
    
    OR4<-(tbl_regression(mod4, exponentiate = TRUE));
    ci4l<-round(OR4$table_body$conf.low[OR4$table_body$label==paste0("CGR")],2);
    ci4h<-round(OR4$table_body$conf.high[OR4$table_body$label==paste0("CGR")],2);  
    est4<-round(OR4$table_body$estimate[OR4$table_body$label==paste0("CGR")],2);
    
   #results
    res4<-c(paste0(x,"=No"),
            est3,
            ci3l,
            ci3h,
            pval1);
    res5<-c(paste0(x,"=Yes"),
            est4,
            ci4l,
            ci4h,NA);
    
    MAT<-matrix(NA,2,5);
    MAT[1,]<-res4;
    MAT[2,]<-res5;
  } 
  return(MAT);
}



customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

IMPT<-function(data_imp,va_dep,va_indep){
  res1<-eval(parse(text=paste0("with(",data_imp,",expr=c('M1'=mean(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[1]]),
                            'SD'=stats::sd(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[1]]) )
                            )"
  )
                            
             )
  )
  res2<-eval(parse(text=paste0("with(",data_imp,",expr=c('M1'=mean(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[2]]),
                            'SD'=stats::sd(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[2]]) )
                            )"
  )
  
  )
  )
  pool_res1<-withPool_MI(res1)
  pool_res2<-withPool_MI(res2)
  
  MAT<-matrix(NA,2,3)
  colnames(MAT)<-c(
    paste0(va_indep),
    paste0(va_dep,"_mean"),
    paste0(va_dep,"_SD")
  )
  
  MAT[1,]<-eval(parse(text=paste0("c(paste(levels(",data_imp,"$data$",va_indep,")[1]),pool_res1[[1]],pool_res1[[2]])"
  )))
  MAT[2,]<-eval(parse(text=paste0("c(paste(levels(",data_imp,"$data$",va_indep,")[2]),pool_res2[[1]],pool_res2[[2]])"
  )))
  MAT[1,1]<-"no transfusion"
  MAT[2,1]<-"transfusion"

  mat<-as.data.frame(MAT)
  mat[,2]<-as.numeric(as.character(mat[,2]))
  mat[,3]<-as.numeric(as.character(mat[,3]))
  print(mat)
}

IMPT_pval<-function(data_imp,va_dep,va_indep){ 
  fit.t.test <- eval(parse(text=paste0("with(data=",data_imp,",exp=lm(",va_dep, "~", va_indep,"))"
)
)
)
t.test.estimates <- pool(fit.t.test)
pval<-summary(t.test.estimates)$p.value[2]
pval1<- ifelse(pval<0.0001, "<0.0001", pval)
print(pval1)
}

IMPT2<-function(data_imp,va_dep,va_indep){
  res1<-eval(parse(text=paste0("with(",data_imp,",expr=c('M1'=mean(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[1]]),
                            'SD'=stats::sd(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[1]]) )
                            )"
  )
  
  )
  )
  res2<-eval(parse(text=paste0("with(",data_imp,",expr=c('M1'=mean(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[2]]),
                            'SD'=stats::sd(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[2]]) )
                            )"
  )
  
  )
  )
  res3<-eval(parse(text=paste0("with(",data_imp,",expr=c('M1'=mean(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[3]]),
                            'SD'=stats::sd(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[3]]) )
                            )"
  )
  
  )
  )
  pool_res1<-withPool_MI(res1)
  pool_res2<-withPool_MI(res2)
  pool_res3<-withPool_MI(res3)
  
  MAT<-matrix(NA,3,3)
  colnames(MAT)<-c(
    paste0(va_indep),
    paste0(va_dep,"_mean"),
    paste0(va_dep,"_SD")
  )
  
  MAT[1,]<-eval(parse(text=paste0("c(paste(levels(",data_imp,"$data$",va_indep,")[1]),pool_res1[[1]],pool_res1[[2]])"
  )))
  MAT[2,]<-eval(parse(text=paste0("c(paste(levels(",data_imp,"$data$",va_indep,")[2]),pool_res2[[1]],pool_res2[[2]])"
  )))
  MAT[3,]<-eval(parse(text=paste0("c(paste(levels(",data_imp,"$data$",va_indep,")[3]),pool_res3[[1]],pool_res3[[2]])"
  )))

  mat<-as.data.frame(MAT)
  mat[,2]<-as.numeric(as.character(mat[,2]))
  mat[,3]<-as.numeric(as.character(mat[,3]))
  mat  
 }

IMPT2_pval<- function(data_imp,va_dep,va_indep){
  pval<-eval(parse(text=paste0("round(mi.anova(mi.res=",data_imp, ",formula='",va_dep, "~", va_indep,"')$anova.table[[1,5]],4)")))
  pval1<- ifelse(pval<0.0001, "<0.0001", pval)
  print(pval1)
  }

miTABLE<-function(data_imp,va_dep,va_indep){ 
  if(eval(parse(text=paste0("is.factor(",data_imp,"$data$",va_dep,")")))){
rescat1Tot <- eval(parse(text=paste0("with(",data_imp,",expr=c('Prop'=prop.table(table(",va_dep,"))[[2]],
                                                              'CIlow'=exactci(table(",va_dep,")[[2]],
                                                             addmargins(table(",va_dep,"))[[3]],
                                                             conf.level = 0.95)[[1]][1],
                                                              'CIup'=exactci(table(",va_dep,")[[2]],
                                                            addmargins(table(",va_dep,"))[[3]],
                                                            conf.level = 0.95)[[1]][2]))")))

res1cat1noT <- eval(parse(text=paste0("with(",data_imp,
                                      ",expr=c(
                                       'Prop'=prop.table(table(",va_indep,",",va_dep,"),1)[[3]],
                                       'CIlow'=exactci(table(",va_indep,",",va_dep,")[[3]],
                                        addmargins(table(",va_indep,",",va_dep,"),2)[[5]],
                                        conf.level = 0.95)[[1]][1],
                                       'CIup'=exactci(table(",va_indep,",",va_dep,")[[3]],
                                       addmargins(table(",va_indep,",",va_dep,"),2)[[5]],
                                       conf.level = 0.95)[[1]][2]
                      ))")))

res1cat1T <- eval(parse(text=paste0("with(",data_imp,
                                    ",expr=c(
                                     'Prop'=prop.table(table(",va_indep,",",va_dep,"),1)[[4]],
                                     'CIlow'=exactci(table(",va_indep,",",va_dep,")[[4]],
                                      addmargins(table(",va_indep,",",va_dep,"),2)[[6]],
                                      conf.level = 0.95)[[1]][1],
                                     'CIup'=exactci(table(",va_indep,",",va_dep,")[[4]],
                                     addmargins(table(",va_indep,",",va_dep,"),2)[[6]],
                                     conf.level = 0.95)[[1]][2]
                    ))")))

Pool_rescat1Tot<-eval(parse(text=paste0("withPool_MI(rescat1Tot)")))
Pool_res1cat1noT<-eval(parse(text=paste0("withPool_MI(res1cat1noT)")))
Pool_res1cat1T<-eval(parse(text=paste0("withPool_MI(res1cat1T)")))

pvalue<-eval(parse(text=paste0("micombine.chisquare(unlist(with(data=",data_imp,",", "exp=chisq.test(",va_indep,",",va_dep,")$statistic)$analyses),df=1)[2]")))
pvalue<-ifelse(pvalue<0.0001, "<0.0001",pvalue)

tab1_cat_tot<-NULL
tab1_cat_tot<-paste(round(Pool_rescat1Tot[["Prop"]],2)*100,
                    "(",
                    round(Pool_rescat1Tot[["CIlow"]],2)*100,
                    "-",
                    round(Pool_rescat1Tot[["CIup"]],2)*100,
                    ")"
)

tab1_cat_noT<-paste(round(Pool_res1cat1noT[["Prop"]],2)*100,
                    "(",
                    round(Pool_res1cat1noT[["CIlow"]],2)*100,
                    "-",
                    round(Pool_res1cat1noT[["CIup"]],2)*100,
                    ")"
)

tab1_cat_T<-paste(round(Pool_res1cat1T[["Prop"]],2)*100,
                  "(",
                  round(Pool_res1cat1T[["CIlow"]],2)*100,
                  "-",
                  round(Pool_res1cat1T[["CIup"]],2)*100,
                  ")"
)
res<-c(va_dep,tab1_cat_tot,tab1_cat_noT,tab1_cat_T,pvalue)
}
else{
  res1Tot <- eval(parse(text=paste0("with(",data_imp, ",expr=c(",'M1',"=mean(",va_dep,"),",'SD',"=stats::sd(",va_dep,")))")))
  
  res1noT <- eval(parse(text=paste0("with(",data_imp, ",expr=c('M1'=mean(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[1]]),
                                   'SD'=stats::sd(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[1]]) ) )")))
  res1T <-   eval(parse(text=paste0("with(",data_imp, ",expr=c('M1'=mean(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[2]]),
                                   SD=stats::sd(",va_dep,"[",va_indep,"==levels(",data_imp,"$data$",va_indep,")[2]]) ) )")))
  
  pool_res1Tot<-eval(parse(text=paste0("withPool_MI(res1Tot)")))
  pool_res1noT<-eval(parse(text=paste0("withPool_MI(res1noT)")))
  pool_res1T<-eval(parse(text=paste0("withPool_MI(res1T)")))
  
  tab_tot1<-paste(round(pool_res1Tot[["M1"]],2),"+/-",round(pool_res1Tot[["SD"]],2))
  tab_noT<-paste(round(pool_res1noT[["M1"]],2),"+/-",round(pool_res1noT[["SD"]],2))
  tab_T<-paste(round(pool_res1T[["M1"]],2),"+/-",round(pool_res1T[["SD"]],2))
  
  fit.t.test1 <- eval(parse(text=paste0("with(data=",data_imp,",", "exp=lm(",va_dep, "~", va_indep,"))")))
  t.test.estimates <- pool(fit.t.test1)
  pval<-round(summary(t.test.estimates)$p.value[2],4)
  pval<-ifelse(pval<0.0001, "<0.0001",pval)
  res<-c(va_dep,tab_tot1,tab_noT,tab_T,pval)
  }
  return(res)
}



INT2<-function(data,var){
  
  dat<-eval(parse(text=paste0("complete(",data,", 'long', inc=TRUE)")))


  dat$x_cat<-as.factor(ifelse(dat[,paste0(var)]>median(dat[,paste0(var)],na.rm=T),paste0(var,">",median(dat[,paste0(var)],na.rm=T)),paste0(var,"<=",median(dat[,paste0(var)],na.rm=T))))
  dat1<-eval(parse(text=paste0("as.mids(dat)")))
  
  mod1<- eval(parse(text=paste0("with(dat1 ,glm(status_1an ~ as.integer(CGR) *x_cat), family = binomial())")))
  OR1<-(tbl_regression(mod1, exponentiate = TRUE));
  
  ci1l<-round(OR1$table_body$conf.low[OR1$table_body$label==paste0("as.integer(CGR)")],2);
  ci1h<-round(OR1$table_body$conf.high[OR1$table_body$label==paste0("as.integer(CGR)")],2);    
  est1<-round(OR1$table_body$estimate[OR1$table_body$label==paste0("as.integer(CGR)")],2);
  pval<-round(OR1$table_body$p.value[[6]],4)
  
  dat$x_cat1<-relevel(dat$x_cat, ref= levels(dat$x_cat)[2]);
  dat2<-eval(parse(text=paste0("as.mids(dat)")))
                                    
  mod2<-eval(parse(text=paste0("with(dat2,glm(status_1an ~ as.integer(CGR) *x_cat1), family = binomial())")))
  
  OR2<-(tbl_regression(mod2, exponentiate = TRUE))
  ci2l<-round(OR2$table_body$conf.low[OR2$table_body$label==paste0("as.integer(CGR)")],2)
  ci2h<-round(OR2$table_body$conf.high[OR2$table_body$label==paste0("as.integer(CGR)")],2)    
  est2<-round(OR2$table_body$estimate[OR2$table_body$label==paste0("as.integer(CGR)")],2)
  
  
  res1<-c(paste0(levels(dat$x_cat)[1]),
          est1,
          ci1l,
          ci1h,
          pval)
  res2<-c(paste0(levels(dat$x_cat1)[1]),
          est2,
          ci2l,
          ci2h,
          "NA")

  MAT<-matrix(NA,2,5)
  MAT[1,]<-res1
  MAT[2,]<-res2
  return(MAT)
  }
 




