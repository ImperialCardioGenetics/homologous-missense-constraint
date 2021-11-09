These are the steps to create HMC scores. The essential scripts and example input files are provided in seperate folders correpsonding to the steps. To run the scripts, the respective file and folder paths would need to be changed.

### Step 1. Seperate each domain alignment into individual files

```python write_stockholm.py $pfam_path/Pfam-A.full.ncbi.gz human_pfam_domain_id_uniq_sort pfam_output_folder error_folder```


### Step 2. From each domain alignment file, we extract the residue: AA, AA position on RefSeq, AA position in the pfam domain alignment (MSA_pos) and pfam_id

sample cmd:
``` python extract_pfam_aa.py a_list_of_pfam_id```
(example input file: step2/human_pfam_domain_id_uniq_sortaa)

### Step 3. Get the synthetic SNVs from the file of Step 2 by merging it with all possible missense variants annotated on human proteome (RefSeq). 

Folder: split_merge_np_missense
```sample cmd: python merge_pfam_domain_np_missense_var.py human_pfam_aa_list```
