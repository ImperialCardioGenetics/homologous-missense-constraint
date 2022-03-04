These are four steps to create HMC scores. The example input and output files by running the steps are provided in the demo (https://github.com/ImperialCardioGenetics/homologous-missense-constraint/tree/main/demo) directory.  


### Step 1. Extract the alignment of Pfam domains annotated in Human

First we download the full alignment file (annotation of protein domain families across species) from Pfam database. From this single big file,  we want to extract the alignment of Pfam domains annotated in Human and save the alignment into an individual file for each domain. This is done by running:

```python write_stockholm.py [Pfam_full_alignment_file] human_pfam_domain_id_list output_folder error_folder```

There are four required arguments: 

[Pfam_full_alignment_file]: Path of full alignment file of all Pfam domains (including human and other species) in stockholm format.
In our study, we used the annotation and full alignment file of all Pfam-A families against NCBI genpept database (Pfam-A.full.ncbi.gz from Pfam version 32.0). It can be downloaded via [Pfam FTP site](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.full.ncbi.gz)). The latest can be found at [current release folder of Pfam FTP site](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release). 

[human_pfam_domain_id_list]: list of human Pfam domain IDs to extract. Example: [step1/human_pfam_domain_id_uniq_sort](https://github.com/ImperialCardioGenetics/homologous-missense-constraint/blob/main/create_HMC/step1/human_pfam_domain_id_uniq_sort).

[output_folder]: path of output folder to store the alignment of pfam domains that we intend to extract

[error_folder]: path of error folder to save the alignment of pfam domains that we exclude


### Step 2. From each domain alignment file, we extract the homologous residue information. 

In the second step, given an input list of pfam domain ids, we extract the homogous residue information for each pfam domain and save in an individual file for each domain. The homologous residue information includes: AA, AA position on RefSeq, AA position in the pfam domain alignment (MSA_pos) and pfam id of the domain. 

This is done by running:

``` python extract_pfam_aa.py [pfam_id_list]```

One required input:

[pfam_id_list]: list of pfam ids to extract homologous residue information. Example: [step2/human_pfam_domain_id_uniq_sortaa](https://github.com/ImperialCardioGenetics/homologous-missense-constraint/blob/main/create_HMC/step2/human_pfam_domain_id_uniq_sortaa)

### Step 3. Get the synthetic SNVs from the file of Step 2 by merging it with all possible missense variants annotated on human transcripts (version: Refseq Select). 

sample cmd: 

```python merge_refseq_select.py list_of_human_pfam_aa_files```

(The input file is a list of files output from step2)

### Step 4: Calculate constraint 

sample cmd:

```Rscript calculate_constraint.R domain_id```
