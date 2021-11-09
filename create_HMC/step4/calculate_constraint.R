

#output file for all domain constraint scores

bug_output="/rds/general/user/xzhang13/home/pfam/analysis_human_pfam/0_generate_constraint/output/calculate_constraint_error.list"

#output dir for domain positions
clean_snv_dir="path_to_output_dir_per_snv"
pos_dir="path_to_output_dir_per_domain_pos"

setwd("~/pfam/analysis_human_pfam")
#correction factor for sites with low sequencing coverage - prepared according to the gnomAD flagship paper
load("/rds/general/user/xzhang13/home/pfam/analysis_human_pfam/0_generate_constraint/coverage_correction_factor.RData")

#mutability from the gnomAD flagship paper - prepared according to the gnomAD flagship paper
load("/rds/general/user/xzhang13/home/pfam/analysis_human_pfam/0_generate_constraint/predicted_prop_AF_0.001.RData")

#AA abreviation hash table
load("/rds/general/user/xzhang13/home/pfam/analysis_human_pfam/0_generate_constraint/AA_abbrev.RData")

args=commandArgs(TRUE)
print(args)
domain_id=args[1]
print(domain_id)


clean_snv_path=paste(clean_snv_dir,domain_id,"_constraint_snp.csv",sep="")
input_path=paste("/rds/general/user/xzhang13/home/pfam/human_pfam_refsel/",domain_id,"_refsel_snp.csv",sep="")
input_file<-read.delim(input_path)

if(nrow(input_file)==0){
  
  #  lines=paste(domain_id,"NA","NA","NA",0,sep="\t")
  #  write(lines,domain_output,append = TRUE)}
}else{
  
  input_file$chr<-sapply(input_file$var,function(x){strsplit(as.character(x),split="_")[[1]][1]})
  input_file$pos<-sapply(input_file$var,function(x){strsplit(as.character(x),split="_")[[1]][2]})
  input_file$nt_change<-sapply(input_file$var,function(x){strsplit(as.character(x),split="_")[[1]][3]})
  input_file$ref<-sapply(input_file$nt_change,function(x){strsplit(x,split="/")[[1]][1]})
  input_file$alt<-sapply(input_file$nt_change,function(x){strsplit(x,split="/")[[1]][2]})
  
  get_column_index<-function(col_name,col_list){
    return(which(col_list==col_name))
  }
  
  
  
  methyl_level<-function(x){
    methy_level = NA
    if(is.na(x)){methy_level=0;return(methy_level)}
    if(x<0.2){methy_level = 0}
    else if(x>0.6){methy_level = 2}
    else {methy_level=1}
    return(methy_level)
  }
  
  get_sub_mutation<-function(x,df){
    col_list=colnames(df)
    met_level=get_column_index("methyl_level",col_list)
    context = get_column_index("three_context",col_list)
    ref = get_column_index("ref",col_list)
    alt = get_column_index("alt",col_list)
    
    sub=NA
    sub=paste(x[context],x[ref],x[alt],sep="_")
    
    return(sub)
  }
  
  #Issue: some variants are in CpG sites but are not methylation mutation. These sites should have zero methlyation level
  #Solution: for all the variants in data frame, if they are not methylation mutation, then set the methylation level to 0
  CpG_mutation<-c("ACG_C_T","CCG_C_T","GCG_C_T","TCG_C_T","CGA_G_A","CGC_G_A","CGG_G_A","CGT_G_A")
  
  check_methylation_level<-function(x,df){
    col_list = colnames(df)
    sub_code<-get_column_index("sub",col_list)
    met_level<-get_column_index("methyl_level",col_list)
    if(!is.element(x[sub_code],CpG_mutation)){
      x[met_level]=0
    }
    return(x[met_level])
  }
  
  #Clean the input such that variants annotated with different NM versions would be reduced 
  #input_file$refseq_NM_id<-sapply(input_file$HGVSc,function(x){strsplit(as.character(x),split="[.]")[[1]][1]})
  #input_file$cDot<-sapply(input_file$HGVSc,function(x){strsplit(as.character(x),split=":")[[1]][2]})
  #input_file$ref_aa<-str_sub(input_file$aa_change,start=3,end=5)
  #input_file<-merge(input_file,AA_abbrev,by.x="ref_aa",by.y="Three",all.x = TRUE)
  #if(sum(is.na(input_file$Abbr))){write(domain_id,bug_output,append=TRUE)}
  
  #input_file=subset(input_file,Abbr==aa)
  select_col=c("pfam_id","domain_len","msa_pos","pfam_AA_seq_id","pfam_AA_seq_len","NP_id","aa_pos","aa","gene_symbol","HGVSc","HGVSp","var","locus","sift","polyphen2","median","methyl_mean","seven_context","three_context","gnomAD_af","chr","pos","nt_change","ref","alt")
  
  uniq_var=unique(input_file[,select_col])
  if(nrow(uniq_var)!=nrow(unique(input_file[,c("pfam_id","msa_pos","var","gene_symbol")]))){write(paste("unique gene msa_pos:",domain_id),bug_output,append=TRUE)}
  clean_file<-subset(uniq_var,duplicated(uniq_var$var)==FALSE)
  if(sum(duplicated(uniq_var$var))>=1){write(paste("unique locus:",domain_id),bug_output,append=TRUE)}
  
  
  
  #Transform the methyl_mean into methyl_level
  clean_file$methyl_level<-sapply(clean_file$methyl_mean,methyl_level)
  
  #Get the substitution code (i.e. context_ref_alt)
  clean_file$sub<-apply(clean_file,1,function(x){get_sub_mutation(x,clean_file)})
  
  #Check the methylation level of the variants, if it's not methylation mutation, the methylation level should be set to 0
  clean_file$methyl_level<-apply(clean_file,1,function(x){check_methylation_level(x,clean_file)})
  
  #Update the substitution code with the right mehtylation level: met_context_ref_alt
  clean_file$mutation<-paste(clean_file$methyl_level,clean_file$sub,sep="_")
  
  #Step 2: merge with substitution rate
  clean_file<-merge(clean_file,predicted_prop,by="mutation",all.x = TRUE)
  
  #input_file$var<-paste(paste(input_file$CHROM,input_file$POS,input_file$REF,sep="_"),input_file$ALT,sep="/")
  #input_file=merge(input_file,id_gnomad_exome_all_chr_mis_AF_0.001,by="var",all.x = TRUE)
  clean_file$observed<-ifelse(clean_file$gnomAD_af!="-",1,0)
  var_file<-subset(clean_file,select = c(var,median,predicted_prop,msa_pos,aa_pos,gnomAD_af,observed))
  
  var_file<-merge(var_file,correction_factor,by.x="median",by.y="coverage",all.x = TRUE)
  
  var_file$adj_exp<-var_file$predicted_prop*var_file$cor_factor_exp
  var_file$exp<-ifelse(!is.na(var_file$adj_exp),var_file$adj_exp,var_file$predicted_prop)
  
  
  #calculate domain constraint and position constraint
  #only consider each genomic variant once 
  
  
  #calculate domain-level constraint and write to the file
  source("/rds/general/user/xzhang13/home/pfam/analysis_human_pfam/domain_score.R")
  domain_score<-domain_constraint(var_file)
  domain_line=paste(domain_id,paste(unlist(domain_score),collapse = "\t"),collapse = "\t")
  clean_file$pfam_ratio<-domain_score$est
  clean_file$pfam_ci_l<-domain_score$ci_l
  clean_file$pfam_ci_u<-domain_score$ci_u
  clean_file$pfam_obs<-domain_score$obs
  clean_file$pfam_exp<-domain_score$exp
  clean_file$pfam_n_var<-domain_score$n_var
  ##### can disable if I don't need all the domain-level constraint scores in one single file 
  #write(domain_line,domain_output,append=TRUE)
  
  #calculate pos constraint and merge with the file
  domain_msa_pos<-sort(unique(var_file$msa_pos))
  msa_pos_map<-data.frame(human_msa_pos=seq(1,length(domain_msa_pos),1),msa_pos=domain_msa_pos)
  var_file<-merge(var_file,msa_pos_map,by="msa_pos",all.x = TRUE)
  pos_score<-sapply(unique(var_file$human_msa_pos),function(x){domain_pos_constraint(var_file,pos=x,k=0)})
  #k1_pos_score<-sapply(unique(var_file$human_msa_pos),function(x){domain_pos_constraint(var_file,pos=x,k=1)})
  #k2_pos_score<-sapply(unique(var_file$human_msa_pos),function(x){domain_pos_constraint(var_file,pos=x,k=2)})
  #k3_pos_score<-sapply(unique(var_file$human_msa_pos),function(x){domain_pos_constraint(var_file,pos=x,k=3)})
  #k5_pos_score<-sapply(unique(var_file$human_msa_pos),function(x){domain_pos_constraint(var_file,pos=x,k=5)})
  #k10_pos_score<-sapply(unique(var_file$human_msa_pos),function(x){domain_pos_constraint(var_file,pos=x,k=10)})
  msa_score<-data.frame(human_msa_pos=unique(var_file$human_msa_pos),residue_ratio=unlist(t(pos_score)[,1]),residue_ci_l=unlist(t(pos_score)[,2]),residue_ci_u=unlist(t(pos_score)[,3]),residue_obs=unlist(t(pos_score)[,4]),residue_exp=unlist(t(pos_score)[,5]),residue_n_var=unlist(t(pos_score)[,6]))
  #                      k1_residue_ci_u=unlist(t(k1_pos_score)[,3]),k1_residue_obs=unlist(t(k1_pos_score)[,4]),k1_residue_exp=unlist(t(k1_pos_score)[,5]),k1_residue_n_var=unlist(t(k1_pos_score)[,6]),
  #                      k2_residue_ci_u=unlist(t(k2_pos_score)[,3]),k2_residue_obs=unlist(t(k2_pos_score)[,4]),k2_residue_exp=unlist(t(k2_pos_score)[,5]),k2_residue_n_var=unlist(t(k2_pos_score)[,6]),
  #                      k3_residue_ci_u=unlist(t(k3_pos_score)[,3]),k3_residue_obs=unlist(t(k3_pos_score)[,4]),k3_residue_exp=unlist(t(k3_pos_score)[,5]),k3_residue_n_var=unlist(t(k3_pos_score)[,6]),
  #                      k5_residue_ci_u=unlist(t(k5_pos_score)[,3]),k5_residue_obs=unlist(t(k5_pos_score)[,4]),k5_residue_exp=unlist(t(k5_pos_score)[,5]),k5_residue_n_var=unlist(t(k5_pos_score)[,6]),
  #                      k10_residue_ci_u=unlist(t(k10_pos_score)[,3]),k10_residue_obs=unlist(t(k10_pos_score)[,4]),k10_residue_exp=unlist(t(k10_pos_score)[,5]),k10_residue_n_var=unlist(t(k10_pos_score)[,6]))
  msa_score<-merge(msa_score,msa_pos_map,by="human_msa_pos",all.x=TRUE)
  clean_file_msa<-merge(clean_file,msa_score,by="msa_pos")
  write.table(clean_file_msa,file=clean_snv_path,sep="\t",col.names = T,row.names = F,quote=F)
  
  domain_pos_output=paste(pos_dir,paste(domain_id,"domain_pos_constraint.txt",sep="_"),sep="/")
  write.table(msa_score,file=domain_pos_output,sep="\t",col.names = T,row.names = F,quote=F)
  
}

