# Arrondi p-valeur
put.pvalue=function(p){
  if (is.na(p)){pc="NA"} else
  if (p<0.0001){pc="<0.0001";} else
    if (0.0001<=p & p<0.001){pc=sprintf("%.4f", p);} else
      if (0.049<=p & p<0.05) {pc=sprintf("%.3f", floor(p*1000)/1000);} else
        if (0.001<=p & p<0.1){pc=sprintf("%.3f", p);} else 
          if (0.1<=p){pc=sprintf("%.2f", p);} 
  return(pc)
}


# Descriptive coninuous variables : Mean ? SD & Median (IQR)

Desc_cont=function(vars,vars.median=NULL,data,round.var=NULL,type=2,
                   verbose=F){
  if (is.null(round.var)){round.var=rep(1,length(vars))};
  sort=NULL;
  N=nrow(data);
  for (k in 1:length(vars)){
    n=sum(!is.na(data[,vars[k]]));
    if (!vars[k]%in%vars.median){
      m=mean(data[,vars[k]],na.rm=T);
      s=sd(data[,vars[k]],na.rm=T);
      temp1=sprintf(paste("%.",round.var[k],"f",sep=""),m)
      temp2=sprintf(paste("%.",round.var[k],"f",sep=""),s)
      temp=c(n,paste(temp1," +/- ",temp2,sep=""));
      sort=rbind(sort,temp);
    } else {
      med=quantile(data[,vars[k]],na.rm=T,type=type,probs=0.5);
      q1=quantile(data[,vars[k]],na.rm=T,type=type,probs=0.25);
      q3=quantile(data[,vars[k]],na.rm=T,type=type,probs=0.75);
      temp1=sprintf(paste("%.",round.var[k],"f",sep=""),med);
      temp2=sprintf(paste("%.",round.var[k],"f",sep=""),q1);
      temp3=sprintf(paste("%.",round.var[k],"f",sep=""),q3);
      temp=c(n,paste(temp1," ","(",temp2," - " ,temp3,")",sep=""));
      sort=rbind(sort,temp);
    }
  }
  rownames(sort)=vars;
  colnames(sort)=c("Available","Mean +/- SD / Median (IQR)");
  if (verbose) {print(sort,quote=F)};
  return(sort);
}


Means2=function(strata,vars,vars.median=NULL,data,
                var.equal=T,vars.nonparam=NULL,
                round.var=NULL,type=2,verbose=F){
  
  if (is.null(round.var)){round.var=rep(1,length(vars));}
  # descritption global;
  N=nrow(data);
  sort=Desc_cont(vars=vars,vars.median=vars.median,
                 data=data,round.var=round.var,
                 type=2);
  
  # description by strata
  data1=data[!is.na(data[,strata]),];
  Nclass=table(data1[,strata]);
  fac=data1[,strata];
  if (!is.factor(fac)){fac=as.factor(fac);}
  cat=levels(fac);
  testNA=rep(F,length(vars))
  for (k in 1:length(cat)){
    temp=Desc_cont(vars=vars,vars.median=vars.median,
                   data=data1[fac==cat[k] & !is.na(fac),],
                   round.var=round.var,type=2);
    sort=cbind(sort,temp);
    testNA=testNA|(temp[,1]==0)
  }
  
  
  # Test
  temp=NULL;
  for (k in 1:length(vars)){
    fml=formula(paste(vars[k],"~as.factor(",strata,")"));
    if (testNA[k]){
      temp=rbind(temp,c("NA","ND"));
    } else {
      if (!vars[k]%in%vars.nonparam){
        test=oneway.test(fml,data=data,var.equal=var.equal)
        pval=test$p.value
        if (var.equal==T){
          temp=rbind(temp,c(put.pvalue(pval),"T-test"));
        } else {
          temp=rbind(temp,c(put.pvalue(pval),"Welch"));
        }
      } else {
        if (length(cat)==2){
          test=wilcox.test(fml,data=data)
          pval=test$p.value
          temp=rbind(temp,c(put.pvalue(pval),"Wilcoxon"));
        } else {
          test=kruskal.test(fml,data=data)
          pval=test$p.value
          temp=rbind(temp,c(put.pvalue(pval),"Kruskal"));
        }
      }
    }
  }
  
  sort=cbind(sort,temp);
  
  # Rename colnames
  col.names=rep(paste("Global"," (N=",N,")",sep=""),2);
  for (j in 1:length(cat)){
    col.names=c(col.names,rep(paste("strata=",cat[j]," (N=",Nclass[j],")",sep=""),2))
  }
  col.names=c(col.names,"p-value","Test");
  colnames(sort)=col.names;
  if (verbose) {print(sort,quote=F)};
  return(sort)
}

# Descriptive qualitative variables / N (percent)

Desc_qual=function(vars,data,desc.all=NULL,round.percent=1,verbose=F){
  N=nrow(data);
  namevar=NULL; 
  sort=NULL;
  for (k in 1:length(vars)){
    tab=table(data[,vars[k]]);
    prop=sprintf(paste("%.",round.percent,"f",sep=""),prop.table(tab)*100);
    if (length(tab)==1){
      namevar=c(namevar,paste(vars[k],"=",names(tab)[1]));
      sort=rbind(sort,c(sum(tab),paste(tab[1]," (",prop[1]," %)",sep="")));
    } 
    if (length(tab)==2 & !(vars[k]%in%desc.all)){
          namevar=c(namevar,paste(vars[k],"=",names(tab)[2]));
          sort=rbind(sort,c(sum(tab),paste(tab[2]," (",prop[2]," %)",sep="")));
    } 
    if (length(tab)>2 | (length(tab)==2 &  (vars[k]%in%desc.all))){
      namevar=c(namevar,paste(vars[k],"=",names(tab)[1]));
      sort=rbind(sort,c(sum(tab),paste(tab[1]," (",prop[1]," %)",sep="")));
      for (j in 2:(length(tab))){
        namevar=c(namevar,paste(vars[k],"=",names(tab)[j]));
        sort=rbind(sort,c("",paste(tab[j]," (",prop[j]," %)",sep="")));
      }
    }
  }
  rownames(sort)=namevar
  colnames(sort)=c("N available","N (%)");
  if (verbose) {print(sort,quote=F)};
  return(sort);
}

# Descriptive + test khi-2

Fisher.test=function(tab){
  tryCatch({
    return(fisher.test(tab));
  },  error = function(e) {
    return(fisher.test(t(tab)));
  })
}

Fisher.chi2=function(tab,correction=F){
  tryCatch({
    return(Fisher.test(tab));
  },  error = function(e) {
    return(chisq.test(tab,correct=correction));
  })
}


Prop2khi2_exact=function(strata,vars,data,desc.all=NULL,
                         chisq.vars=NULL,exact.vars=NULL,
                         round.percent=1,verbose=F){
  data=data[!is.na(data[,strata]),];
  sort=varname=NULL;
  N=length(which(is.na(data[,strata])==F))
  Nclass=table(data[,strata]); 
  categ=names(Nclass)
  for (k in 1:length(vars)){
    tab=table(data[,vars[k]],data[,strata])
    prop=matrix(sprintf(paste("%.",round.percent,"f",sep=""),prop.table(tab,2)*100),nrow=nrow(tab),ncol=ncol(tab),byrow=F);
    class=paste(vars[k],"=",rownames(tab));
    tab.expected=outer(rowSums(tab), colSums(tab), "*")/sum(tab)
    testDone=F;
    
    if (testDone==F & vars[k]%in%chisq.vars){      
      if (nrow(tab)==1){
        p.chisq=NA;
        pc="NA";
        type.test="ND";
      } else  {
        if (min(tab.expected)<5){correction=T;} else {correction=F}
        p.chisq=chisq.test(tab,correct=correction)$p.value;
        if (!(is.na(p.chisq))) {
          pc=put.pvalue(p.chisq);
          type.test="Chi-square";
          if (correction==T){type.test="Chi-square (with continuity correction)";}
        } else {
          pc="NA";
          type.test="ND";
        }
      }
      testDone=T;
    }
    
    if (testDone==F & vars[k]%in%exact.vars){      
      if (nrow(tab)==1){
        p.chisq=NA;
        pc="NA";
        type.test="ND";
      } else  {
        if (min(tab.expected)<5 & max(dim(tab))<=2){correction=T;} else {correction=F}
        temp=Fisher.chi2(tab,correction=correction)
        p.exact=temp$p.value;
        pc=put.pvalue(p.exact);
        if (temp$method=="Fisher's Exact Test for Count Data") {type.test="Fisher's exact";} else {type.test="Chi-square";}
      }
      testDone=T;
    }
    
    if (testDone==F){
      if (nrow(tab)==1){
        p.chisq=NA;
        pc="NA";
        type.test="ND";
      } else  {
        if (min(tab.expected)>=5){
          p.chisq=chisq.test(tab,correct=F)$p.value;
          pc=put.pvalue(p.chisq);
          type.test="Chi-square";
        } else {
          p.exact=fisher.test(tab)$p.value;
          pc=put.pvalue(p.exact);
          type.test="Fisher's exact";
        }
      }
      testDone=T;
    }
    
    if (nrow(tab)==1){Navailable=tab} else {Navailable=apply(tab,2,sum);}
    
    temp=NULL;
    
    if (length(class)==1){
      for (l in 1:length(categ)){
        temp=c(temp,c(Navailable[l],paste(tab[l]," (",prop[l]," %)",sep="")));
      }
      sort=rbind(sort,c(class[1],temp,pc,type.test));
      varname=c(varname,class[1]);
    } 
    
    if (length(class)==2 & !(vars[k]%in%desc.all)){
      for (l in 1:length(categ)){
        temp=c(temp,c(Navailable[l],paste(tab[2,l]," (",prop[2,l]," %)",sep="")));
      }
      sort=rbind(sort,c(class[2],temp,pc,type.test));
      varname=c(varname,class[2]);
    } 
    
    if ((length(class)>2) | (length(class)==2 &  (vars[k]%in%desc.all))){
      
      for (j in 1:length(class)){
        
        if (j>=2){
          temp=NULL;
          for (l in 1:length(categ)){
            temp=c(temp,c("",paste(tab[j,l]," (",prop[j,l]," %)",sep="")))
          }
          sort=rbind(sort,c(class[j],temp,"",""));
        } else {
          temp=NULL;
          for (l in 1:length(categ)){
            temp=c(temp,c(Navailable[l],paste(tab[j,l]," (",prop[j,l]," %)",sep="")))
          }
          sort=rbind(sort,c(class[j],temp,pc,type.test));
        }
      }
      varname=c(varname,class);
    }
  }
  
  col.names=NULL;
  for (j in 1:length(categ)){
    col.names=c(col.names,paste("(N=",Nclass[j],")",sep=""),paste(strata,"=",categ[j],sep=""))
  }
  colnames(sort)=c("Categories",col.names,"p-value","Test");
  if (verbose) {print(sort,quote=F);}
  return(sort);
}


# Table 1
Table1=function(strata,vars,vars.cont=NULL,vars.qual=NULL,desc.all=NULL,
                round.var.cont=NULL,vars.median=NULL,
                var.equal=T,vars.nonparam=NULL,
                chisq.vars=NULL,exact.vars=NULL,
                round.percent=1,type=2,verbose=F,
                data,put.effectif=T){
  sort=NULL;varnames=NULL;
  if (is.null(round.var.cont)){round.var.cont=rep(1,length(vars.cont))};
  j=1;
  for (k in 1:length(vars)){
    if (vars[k]%in%vars.cont){
      temp=Means2(strata,vars=vars[k],vars.median=vars.median,data=data,
                  var.equal=var.equal,vars.nonparam=vars.nonparam,
                  round.var=round.var.cont[j],type=2,verbose=F);
      varnames=c(varnames,vars[k]);
      sort=rbind(sort,temp);
      j=j+1;
    }
    
    if (vars[k]%in%vars.qual){
      temp1=Desc_qual(vars=vars[k],data=data,desc.all=desc.all,
                      round.percent=round.percent);
      temp2=Prop2khi2_exact(strata,vars=vars[k],data=data,desc.all=desc.all,
                            chisq.vars=chisq.vars,exact.vars=exact.vars,
                            round.percent=round.percent);
      temp3=cbind(temp1,temp2)[,-3];
      varnames=c(varnames,temp2[,1]);
      sort=rbind(sort,temp3);
    }
    
  }
  rownames(sort)=varnames;
  tab=table(data[,strata]);
  prop=sprintf("%.1f",prop.table(tab)*100);
  temp1=NULL;
  for (k in 1:length(tab)){
    temp1=c(temp1,paste(strata,"=",names(tab)[k],sep=""),
            paste(strata,"=",names(tab)[k],sep=""));
  }
  colnames(sort)=c("Global",paste("Global", " (N=",sum(tab),")",sep=""),temp1,"p-value","Test");
  if (verbose) {print(sort,quote=F);}
  
  if (put.effectif==F){
    sort=sort[,-seq(3,ncol(sort)-3,2)];
  }
    
  return(sort)
}


Desc_global=function(vars,varcont,round.var,vars.desc.all,data){
  tab1=NULL;  
  n=1;
  for (k in 1:length(vars)){
    if (vars[k]%in%varcont){
      temp1=Desc_cont(vars=vars[k],vars.median=NULL,data=data,round.var=round.var[n],type=2,
                      verbose=F);
      temp2=Desc_cont(vars=vars[k],vars.median=vars[k],data=data,round.var=round.var[n],type=2,
                      verbose=F);
      tab1=rbind(tab1,cbind(temp1,temp2[,2]));
      n=n+1;
    } else {
      temp1=Desc_qual(vars=vars[k],data=data,desc.all=vars.desc.all,round.percent=1,verbose=F)
      tab1=rbind(tab1,cbind(temp1,""))
    }
  }
  colnames(tab1)=c("Available","Mean ? SD / n (%)","Median (IQR)")
  print(tab1,quote=F)
  return(tab1)
}


