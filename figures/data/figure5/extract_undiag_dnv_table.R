#Rscript to create the table of de novo variants observed in the undiagnostic probands

original <- read_excel("Desktop/OneDrive/OneDrive - Imperial College London/domain_constraint/2020_Kaplanis/41586_2020_2832_MOESM4_ESM_S2.xlsx")
all_dnv_update_31KDD <- read.delim("~/Desktop/OneDrive/OneDrive - Imperial College London/domain_constraint/2020_Kaplanis/20211201/input/update_denovowest_20211201/all_dnv_update_31KDD.txt", stringsAsFactors=FALSE)

diagnostic_gene<-original[which(original$diagnostic_category=="consensus"),"symbol"]$symbol

cases_withdiagnosedgene<-subset(all_dnv_update_31KDD,is.element(symbol,diagnostic_gene))

syn_var<-c("synonymous_variant","coding_sequence_variant")
cases_withdiagnosedvar<-subset(all_dnv_update_31KDD,is.element(symbol,diagnostic_gene)&!is.element(cq,syn_var))
diagnostic_cases=unique(cases_withdiagnosedvar$id)

all_dnv_update_31KDD$diagnosed_cases<-ifelse(is.element(all_dnv_update_31KDD$id,diagnostic_cases),TRUE,FALSE)
table(all_dnv_update_31KDD$diagnosed_cases)

undiagnosed_cases<-subset(all_dnv_update_31KDD,!is.element(id,diagnostic_cases))

##sanity check
#number of cases with DNV in the full cohort
length(unique(all_dnv_update_31KDD$id))
## number of undiagnosed cases
length(unique(undiagnosed_cases$id))
## number of diagnosed cases - 6770
length(unique(diagnostic_cases))


write.table(undiagnosed_cases,file="~/Desktop/OneDrive/OneDrive - Imperial College London/domain_constraint/2020_Kaplanis/20211201/input/update_denovowest_20211201/undiag_dnv_24KDD_20211205.txt",sep="\t",quote=F,col.names = T,row.names=F)
