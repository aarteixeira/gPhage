# gPhage NGS analysis
Analyze genomic phage display NGS data.

Teixeira AAR, et al. A refined genome phage display methodology delineates the human antibody response in patients with Chagas disease. iScience 2021; 24:102540. 

## dna_processing.py
Process fastq files to extract insert sequences and compare to reference genomes and proteomes.

Customize variables to adjust to you pipeline:
| Variable | Description |
|:---:|:---:|
|BLASTN_EXE | path/command to BlastN executable|
|BLASTN_DB | path to BlastN database |
|BLASTP_EXE | 	path/command to BlastP executable |
|BLASTP_DB |   	path to BlastP database |
|dict_sample_fastq | Sample names and path to respective FastQ to be processes. To add a sample, add a new line with the following code: dict_sample_fastq["Sample Name"] = "path_to_fastq.fastq" |
|cutoff_dna_identity | Minimum relative identity to reference genome for a DNA insert to be considered a hit (default is 0.9) |
|cutoff_peptide_similarity | Minimum relative similarity to reference proteome for a peptide to be considered a hit (default is 0.6) |
|flank_upstream  | sequence found upstream of the DNA insert (default set for pG8SAET sequence) |
|flank_downstream | sequence found downstream of the DNA insert (default set for pG8SAET sequence) |


The script will generate the following output:
|Output File|	Description|
|:---:|:---:|
|dna_dataframe.xlsx|  	table containing all DNA inserts identified in the all samples, respective frequencies in each sample, including the frequency in each sample and BlastN results against reference genomes.|
|dna_above_cutoff.xlsx| 	same as dna_dataframe.xlsx, but containing only the inserts that are above the BlastN identity cutoff (match reference genomes)|
|dna_below_cutoff.xlsx|   	same as dna_dataframe.xlsx, but containing only the inserts that are below the BlastN identity cutoff (do not match reference genomes) |
|peptide_dataframe.xlsx|	table containing all peptides identified in the all samples peptide, including the frequency in each sample and BlastP results against reference proteomes.|
|norm_peptide_df.xlsx|	same as peptide_dataframe.xlsx, but with counts normalized by total DNA (relative frequency)|
|peptide_above_cutoff.xlsx|	same as peptide_dataframe.xlsx, but containing only the peptides that are above the BlastP similarity cutoff (match reference proteomes)|
|peptide_below_cutoff.xlsx|	same as peptide_dataframe.xlsx, but containing only the peptides that are below the BlastP similarity cutoff (do not match reference proteomes)|
|statistics.xlsx|	general processing statistics|
|peptides.txt|	A list of all peptides identified by the script. This file can used as input for the clustering.py script|


### Packages required:
* Pandas (v1.3.0)
* Biopython (v1.79)
* Numpy (v1.21.0)


## clustering.py
Cluster peptide sequences to identify epitopes.

Customize variables to adjust to you pipeline:
|Variable|	Description|
|:---:|:---:|
|PEPTIDES_TXT|	Path to the file containing the peptides to be clustered (text file containing one peptide sequence per line)|
|MAFFT|	path/command to MAFFT executable|
|HMMBUILD|	path/command to hmmbuild executable|
|HMMEMIT|	path/command to hmmemit executable|
|K|	K-mer length for peptide search. Only peptides that share a K-mer are compared. Default is 4 amino acids.|
|min_len|	Minimum length of the peptides used for clustering. Peptides with length below this threshold are excluded. Default is 6 amino acids.|
|score_cutoff|	Minimum score (% partial identity) for peptides to be clustered together. Default is 80.|

The script will generate the following output:
|Output File|	Description|
|:---:|:---:|
|cluster_consensus.xlsx|	table containing cluster IDs, number of different peptides forming the clusters, and consensus sequence.|
|peptide_cluster.xlsx| 	table containing all peptide sequences and corresponding cluster.|
|align folder|	Folder containing the multiple sequence alignments, consensus sequences, and HMM profiles for each cluster.|

### Packages required:
* Pandas (v1.3.0)
* Biopython (v1.79)
* Fuzzywuzzy (v0.18.0)


For questions: andrearteixeira@gmail.com