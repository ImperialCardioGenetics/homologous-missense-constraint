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

### Step 3. Annotate all the possible missense variants at the homologous residues 

In the third step, we get all the possible missense variants at the homologous residues by:
(1)synthesize all possible single nucleotide variants (SNVs) at the homologous residues
(2)keep the SNVs that are annotated as missense variants on human transcripts (in our study, we used Refseq Select as reference)

This is done by running:

```python merge_refseq_select.py [all_human_mis_snv] [human_pfam_aa_file_list]```

Two required inputs:

[all_human_mis_snv]: a flat file with all possible missense variants annotated at human transcripts. Example: [step3/example_all_pos_rare_mis_refseq.txt](https://github.com/ImperialCardioGenetics/homologous-missense-constraint/blob/main/create_HMC/step3/example_all_pos_rare_mis_refseq.txt)

[human_pfam_aa_file_list]: a list of homologous residue information files (output file from Step 2). Example: [step3/sub_human_pfam_aa_list_aa](https://github.com/ImperialCardioGenetics/homologous-missense-constraint/blob/main/create_HMC/step3/sub_human_pfam_aa_list_aa)

### Step 4: Calculate genetic constraint for each homologous residues 

This is done by running: 

(??change the R script: input arguments and output path)

```Rscript calculate_constraint.R [input_domain_homologous_mis] [coverage_correction_facotor] [mutability_table] [AA_abreviation_table]```

Required inputs: 

[input_domain_homologous_mis]: a file with all possible missense variants at the homologous residues in a protein domain family (output file from Step 3)

The other required inputs are assistant files to calculate constraint, we have provided the ones used in our study:

[coverage_correction_factor]: An R data object saving the factors used to adjust the probability of neutral substitions for genomic sites with low sequencing-coverage (More information can be found at the [Zhang.X, et.al medRxiv](https://www.medrxiv.org/content/10.1101/2022.02.16.22271023v1) Method section - Developing a selection-neutral, sequence-context mutational model). 

[mutability_table]: An R data object saving the probabilities of neutral substitutions in gnomAD 125K exomes for each tri-nucleotide sequence context and methylation level. This is calculated by calibration from baseline mutation rate to probabilities of synonymous variants in gnomAD 125K exomes.(More information can be found at the [Zhang.X, et.al medRxiv](https://www.medrxiv.org/content/10.1101/2022.02.16.22271023v1) Method section - Developing a selection-neutral, sequence-context mutational model)

(?? more specific)
[AA_abreviation_table]: An R data object saving the table of amino acid abreviation  
