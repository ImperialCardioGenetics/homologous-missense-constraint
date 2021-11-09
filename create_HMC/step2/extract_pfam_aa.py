import sys
import os.path

def extract_seq(pfam_id,dir_path):
    
    pfam_file=open(os.path.join(dir_path,pfam_id+".txt"),"r")
    line = pfam_file.readline()
    ids={}
    seqs={}
    passed_end_alignment = False
    
    while line and not passed_end_alignment: #find the line where the query number would be there
            
        if line[:2] == "//": #the end of the alignment and leave the loop 
            passed_end_alignment = True
            break;
        elif line == "": # blank line, ignore
            pass
        elif line == "\n":
            pass
        elif line[0] != "#":
            # Sequence
            # Format: "<seqname> <sequence>"
            assert not passed_end_alignment
            parts = [x.strip() for x in line.split(" ", 1)]
            if len(parts) != 2:
                print(parts)
                # This might be someone attempting to store a zero length sequence?
                raise ValueError(
                    "Could not split line into identifier "
                    "and sequence:\n" + line)
            seq_id, seq = parts
            if seq_id not in ids:
                ids[seq_id] = True
                seqs.setdefault(seq_id, '')
                seqs[seq_id] += seq.replace(".", "-") 
        line=pfam_file.readline()
        
    assert len(seqs) <= len(ids)
    
    output_path=dir_path+"_cmd"
    output=open(os.path.join(output_path,pfam_id+"_AA.txt"),"w")
    output.write('pfam_id'+"\t"+'domain_len\t'+'msa_pos'+'\t'+'pfam_AA_seq_id\t'+'pfam_AA_seq_len\t'+'NP_id'+"\t"+'aa_pos'+"\t"+'aa'+"\n")

    for seq_id in ids:
            #only get the human protein sequences
            #get the start and end pos of the protein sequence
        refseq_id = seq_id.split("/")[0]
        NP = refseq_id.startswith('NP')
        if NP:
            start = int(seq_id.split("/")[1].split("-")[0])
            end = int(seq_id.split("/")[1].split("-")[1])
            seq = seqs[seq_id]
            aa_pos = start
            AA_len = sum(x.isalpha() for x in seq)
            for i in range(len(seq)):
                msa = i + 1
                    # 1-coordinate system    
                if seq[i].isalpha():
                	output.write(pfam_id+"\t"+str(len(seq))+"\t"+str(msa)+"\t"+seq_id+"\t"+str(AA_len)+"\t"+refseq_id+"\t"+str(aa_pos)+"\t"+seq[i]+"\n")
                	aa_pos=aa_pos+1
                
            assert aa_pos - 1 == end
            assert msa <= len(seq)
            
    output.close()
    return 1



human_pfam_file_path=sys.argv[1]
pfam_list = [line.rstrip() for line in open(human_pfam_file_path)]

dir_path = "$pfam_output_folder"
for pfam_id in pfam_list:
    print(pfam_id)
    extract_seq(pfam_id,dir_path)
