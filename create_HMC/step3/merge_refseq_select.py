import pandas as pd
import os, fnmatch
import sys

mis_file=pd.read_csv("/rds/general/user/xzhang13/projects/lms-ware-analysis/live/xzhang/all_pos_SNV/refseq_select_seq_context/refseq_select_rare_mis_cont.txt",sep="\t",low_memory=True)
human_pfam_file_path = sys.argv[1]

pfam_list = [line.rstrip() for line in open(human_pfam_file_path)]
for file in pfam_list:
    pfam_file = pd.read_csv(file,sep="\t")
    if(pfam_file.shape[0]):
        pfam_id=pfam_file['pfam_id'][0]
        new_file=pfam_file.merge(mis_file,left_on=['NP_id','aa_pos'],right_on=['NP_id','AA_pos'])
        if(new_file.shape[0]):
            new_file.to_csv('/rds/general/user/xzhang13/home/pfam/human_pfam_refsel/'+pfam_id+"_refsel_snp.csv", sep="\t",index=False)
            print(pfam_id+"\n")
