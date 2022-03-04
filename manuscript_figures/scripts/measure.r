ppv<-function(th,df,score){
  df[,"label"]<-as.numeric(as.character(df[,"label"]))
  pred_positive=df[which(df[,score]<=th),]
  return(sum(pred_positive$label)/nrow(pred_positive))
}

tpr<-function(th,df,score){
  df[,"label"]<-as.numeric(as.character(df[,"label"]))
  positive<-subset(df,label==1)
  pred_label<-ifelse(positive[,score]<=th,1,0)
  return(sum(pred_label)/nrow(positive))
}

tpr_bin<-function(th1,th2,df,score){
  df[,"label"]<-as.numeric(as.character(df[,"label"]))
  positive<-subset(df,label==1)
  pred_label<-ifelse(positive[,score]<th2&positive[,score]>=th1,1,0)
  return(sum(pred_label)/nrow(positive))
}

tnr<-function(th,df,score){
  df[,"label"]<-as.numeric(as.character(df[,"label"]))
  negative<-subset(df,label==0)
  pred_label<-ifelse(negative[,score]<th,1,0)
  return(1-(sum(pred_label)/nrow(negative)))
}


#for the benign variants 
fpr<-function(th,df,score){
  df[,"label"]<-as.numeric(as.character(df[,"label"]))
  pred_positive=df[which(df[,score]<=th),]
  false<-subset(pred_positive,label==0)
  negative=df[which(df$label==0),]
  return(nrow(false)/nrow(negative))
}

fpr_bin<-function(th1,th2,df,score){
  df[,"label"]<-as.numeric(as.character(df[,"label"]))
  pred_positive=df[df[,score]<th2&df[,score]>th1,]
  false<-subset(pred_positive,label==0)
  negative=df[which(df$label==0),]
  return(nrow(false)/nrow(negative))
}


ppv_scale<-function(th,scale,df,score){
  df[,"label"]<-as.numeric(as.character(df[,"label"]))
  df=df[which(df[,score]>th),]
  df=df[which(df[,score]<(th+scale)),]
  pred_label<-df[which(df[,score]<th),]
  return(sum(pred_label$label)/nrow(pred_label))
}

tpr_scale<-function(th,scale,df,score){
  df[,"label"]<-as.numeric(as.character(df[,"label"]))
  df=df[which(df[,score]>th),]
  df=df[which(df[,score]<(th+scale)),]
  positive<-subset(df,label==1)
  pred_label<-ifelse(positive[,score]<th,1,0)
  return(sum(pred_label)/nrow(positive))
}



mcc<-function(th,scale,df,score){
  df=df[which(df[,score]>th),]
  df=df[which(df[,score]<(th+scale)),]
  true<-subset(df,label==1)
  false<-subset(df,label==0)
  tp<-1.0*sum(true[,score]<th)
  fp<-1.0*nrow(true)-tp
  tn<-1.0*sum(false[,score]>th)
  fn<-1.0*nrow(false)-tn
  mcc=(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  return(mcc)
}

scale_confusion_matrix<-function(th,scale,df,score){
  df[,'label']<-as.numeric(as.character(df[,"label"]))
  df=df[which(df[,score]>th),]
  df=df[which(df[,score]<(th+scale)),]
  true<-subset(df,label==1)
  false<-subset(df,label==0)
  tp<-1.0*sum(true[,score]<th)
  fp<-1.0*nrow(true)-tp
  tn<-1.0*sum(false[,score]>th)
  fn<-1.0*nrow(false)-tn
  mat<-matrix(c(tp,fp,fn,tn),2,2,byrow=TRUE)
  
  return(mat)
}

confusion_matrix<-function(th,df,score){
  df[,'label']<-as.numeric(as.character(df[,"label"]))
  true<-subset(df,label==1)
  false<-subset(df,label==0)
  tp<-1.0*sum(true[,score]<th)
  fp<-1.0*nrow(true)-tp
  tn<-1.0*sum(false[,score]>th)
  fn<-1.0*nrow(false)-tn
  if(fp==0|tp==0|fn==0|tn==0){fp=fp+0.5;tn=tn+0.5;tp=tp+0.5;fn=fn+0.5}
  mat<-matrix(c(tp,fp,fn,tn),2,2,byrow=TRUE)
  
  return(mat)
}

odds_ratio<-function(mat){
  or=(mat[1]/mat[3])/(mat[2]/mat[4])
  se_log_or<-sqrt(1/mat[1]+1/mat[2]+1/mat[3]+1/mat[4])
  ci_u=exp(log(or)+1.96*se_log_or)
  ci_l=exp(log(or)-1.96*se_log_or)
  return(c(or,ci_l,ci_u))
  
  }
  
confusion_matrix_range<-function(df,score,th_l,th_u){
  ## this is depleted as it's problematic to include smaller constraint as TN
  ############
  ###                 Pathogenic            Benign
  ### in the bin (low~high)      
  ### out of the bin (<=low|>=high)
  
  df[,'label']<-as.numeric(as.character(df[,"label"]))
  true<-subset(df,label==1)
  false<-subset(df,label==0)
  tp<-1.0*sum(true[,score]<=th_u&true[,score]>th_l)
  #tp<-1.0*sum(true[,score]<th)
  fn<-1.0*nrow(true)-tp
  tn<-1.0*sum(false[,score]>=th_u | false[,score]<th_l)
  fp<-1.0*nrow(false)-tn
  if(fn==0){fn=0.5;tp=tp+0.5;tn=tn+0.5;fp=fp+0.5}
  mat<-matrix(c(tp,fp,fn,tn),2,2,byrow=TRUE)
  
  return(mat)
}
library(fmsb)
enrichment_range<-function(df,score,th_l,th_u){
  df["label"]<-as.numeric(as.character(df[,"label"]))
  true<-subset(df,label==1)
  false<-subset(df,label==0)
  tp<-1.0*sum(true[,score]<=th_u&true[,score]>th_l)
  fp<-1.0*sum(false[,score]<=th_u&false[,score]>th_l)
  fn<-nrow(true)-tp
  tn<-nrow(false)-fp
  
  outcome<-riskratio(tp,fp,tp+fn,fp+tn)
  
  
  return(c(outcome$estimate,outcome$conf.int))
}

## so that p.v will print precise decimal values
log10pvalue_riskratio<-function(df,score,th_l,th_u,lower.tail=FALSE){
  df["label"]<-as.numeric(as.character(df[,"label"]))
  true<-subset(df,label==1)
  false<-subset(df,label==0)
  tp<-1.0*sum(true[,score]<=th_u&true[,score]>th_l)
  fp<-1.0*sum(false[,score]<=th_u&false[,score]>th_l)
  fn<-nrow(true)-tp
  tn<-nrow(false)-fp
  X=tp
  Y=fp
  m1=tp+fn
  m2=fp+tn
  ESTIMATE <- (X/m1)/(Y/m2)
  n1 <- X + Y
  Total <- m1 + m2
  n2 <- Total - n1
  p.v <-  log10(exp(1))*(log(2)+pnorm(abs((X - n1 * m1/Total)/sqrt(n1 * n2 * m1 * m2/Total/Total/(Total - 1))),lower.tail = lower.tail,log.p=TRUE))
  return(p.v)
}




