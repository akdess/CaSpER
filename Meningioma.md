---
title: "Bulk Meningioma DataSet"
output: github_document
---

CaSpER is an algorithm for identification, visualization and integrative analysis of CNV events in multiscale resolution using single-cell or bulk RNA sequencing data.

It takes as input: 
- Mapped RNA-Seq reads (Bam files)
- Normalized expression matrix

and outputs: 
- CNv events in mutliscale resolution
- Mutually exclusive and co-occurent CNV events


Installation
----------

CaSpER is developed under R version 3.4.3.  

You can install CaSpER R package using the following R commands:

``` r
require(devtools)
install_github("akdess/CaSpER")

```

For extracting B-allele freqeuncies from RNA-Seq bam files download BAFExtract c++ source or binary code from here. The easiest way is to download the binary file.


After downloading  source code type the following: 
``` r
cd loh
make clean
make

```
The executable is located under directory /bin. 

The  docker version of BAF extract is available through XXX. 
 

Usage
----------

Step 1. Merge single cell RNA-Seq bam files (required for single-cell studies)

```{bash} 
	 /home1/05227/akdes/bin/bamtools merge -list <bam_file_names> -out <sample_name>_merged.bam 
	 samtools index <sample_name>_merged.bam
```

Step 2. Extract BAF values from RNA-Seq bam files
	
```{bash}
	samtools view <bam_file> | ./loh -generate_compressed_pileup_per_SAM stdin /home1/05227/akdes/work/hg19.list <sample_dir> 50 0; ./mapped_read_tools -get_SNVs_per_pileup  <genome_list> <output_dir> <genome_fasta_dir> 20 4 0.1 <output_file>
```

To generate <genome_fasta_dir>: 
	
hg38: 
```{bash}
		wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz; tar -xvzf hg38.chromFa.tar.gz
		cd <genome_fasta_dir>; genome_sequence_tools -binarize_fasta_per_dir . . fa
```
hg19: 
```{bash}
		wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz; tar -xvzf hg38.chromFa.tar.gz
		cd <genome_fasta_dir>; genome_sequence_tools -binarize_fasta_per_dir . . fa
```
	
To generate genome_list file: 
download http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes

hg38:
```{bash}
	fetchChromSizes hg38 > hg38.list
```
	
hg19: 
```{bash}
	fetchChromSizes hg19 > hg19.list
```


Step 3. Run CaSpER R package (See tutorials below)


Tutorials
----------

1. Yale meningioma Bulk RNA-Seq dataset [another](Meningioma.md)
2. TCGA-GBM Bulk RNA-Seq dataset
3. GBM Single-cell RNA-Seq dataset