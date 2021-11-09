import gzip
import sys
from collections import OrderedDict
import os.path

def write_stockholm(file_handle,domain_list,human_path,other_path):
    line = file_handle.readline()
    i=0
    domain_find = {id:False for id in domain_list}
    
    
    while line: #find the line where the query number would be there
        
            
        if line[:7] == "#=GF AC":
            AC = line[7:].strip().split(".",1)[0]
            line_set=[]
            end_alignment = False
            
            while(line and not end_alignment): 
                if line[:2] == "//":
                    end_alignment = True
                line_set.append(line)
                line=file_handle.readline()
                
            if(AC):
                sys.stdout.write("human_pfam:"+AC+'\n')
                if(AC in domain_list and domain_find[AC]==False):
                    path = os.path.join(human_path,AC+".txt")
                    output = open(path,"w")
                    for line in line_set:
                        output.write(line)
                    output.close()
                else:
                    path = os.path.join(other_path,AC+".txt")
                    output = open(path,"w")
                    for line in line_set:
                        output.write(line+"\n")
                    output.close()
            else:
                sys.stdout.write("other:"+str(i)+'\n')
                path = os.path.join(other_path,"check_"+str(i)+".txt")
                output = open(path,"w")
                for line in line_set:
                    output.write(line)
                output.close()
                i=i+1
        else: 
            line=file_handle.readline()
        
                

if __name__ == '__main__':
    pfam_file = sys.argv[1] #by default assume it's in gzip format
    domain_file = sys.argv[2]
    human_path = sys.argv[3] #output path
    other_path = sys.argv[4]
    pfam = gzip.open(pfam_file,"rt",encoding='latin-1')
    domain_list = [line.rstrip() for line in open(domain_file)]
    
    pfam = write_stockholm(pfam, domain_list,human_path,other_path)
