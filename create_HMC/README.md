These are the steps to create HMC scores. The essential scripts and example input files are provided in seperate folders correpsonding to the steps. To run the scripts, the respective file and folder paths would need to be changed.

### Step 1. Seperate each domain alignment into individual files

```python write_stockholm.py $pfam_path/Pfam-A.full.ncbi.gz human_pfam_domain_id_uniq_sort pfam_output_folder error_folder```


### Step 2. From each domain alignment file, we extract the residue: AA, AA position on RefSeq, AA position in the pfam domain alignment (MSA_pos) and pfam_id

sample cmd:

``` python extract_pfam_aa.py a_list_of_pfam_id```

(example input file: step2/human_pfam_domain_id_uniq_sortaa)

### Step 3. Get the synthetic SNVs from the file of Step 2 by merging it with all possible missense variants annotated on human transcripts (version: Refseq Select). 

sample cmd: 

```python merge_refseq_select.py list_of_human_pfam_aa_files```

(The input file is a list of files output from step2)

### Step 4: Calculate constraint 

sample cmd:

```Rscript calculate_constraint.R domain_id```
