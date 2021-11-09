domain_constraint<-function(df,alpha=0.05){
  obs<-sum(df$observed,na.rm=TRUE)
  exp<-sum(df$exp,na.rm = TRUE)
  constraint_ci<-constraint_ci(obs,exp,alpha)
  ratio<-obs/exp
  ci_l<-constraint_ci[1]
  ci_u<-constraint_ci[2]
  return(list(est=ratio,ci_l=ci_l,ci_u=ci_u,obs=obs,exp=exp,n_var=nrow(df)))
}

domain_pos_constraint<-function(df,pos,k=0,alpha=0.05){
  df_pos<-subset(df,human_msa_pos<=pos+k&human_msa_pos>=pos-k)
  return(domain_constraint(df_pos,alpha))
}


#using bayesian oe_constraint estimate in gnomAD paper 
#credit: https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/constraint_basics.py 
#oe_constraint_interval

constraint_ci<-function(obs,exp,alpha=0.05){
  l=seq(1,2000,1)/1000
  range_p<-sapply(l,function(x){dpois(obs,exp*x)})
  
  cum_sum_arr=c()
  cum_sum=0
  for(i in 1:length(range_p)){
    cum_sum=cum_sum+range_p[i]  
    cum_sum_arr=append(cum_sum_arr,cum_sum)
  }
  
  
  cumulative_p<-cum_sum_arr
  max_cumulative_dpois=max(cumulative_p)
  norm_dpois=cumulative_p/max_cumulative_dpois
  
  lower_idx=ifelse(obs>0,max(which(norm_dpois<alpha)),0)
  upper_idx=min(which(norm_dpois>1-alpha))
  oe_ci_l=lower_idx/1000
  oe_ci_u=upper_idx/1000
  return(c(oe_ci_l,oe_ci_u))
 }

