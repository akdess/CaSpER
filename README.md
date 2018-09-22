CaSpER
---

CaSpER is an algorithm for identification, visualization and integrative analysis of CNV events in multiscale resolution using single-cell or bulk RNA sequencing data.

It takes as input: 

- Mapped RNA-Seq reads (Bam files)

- Normalized expression matrix

and outputs: 

- CNV events in multiscale resolution

- Mutually exclusive and co-occurent CNV events


Installation
----------

CaSpER is developed under R version 3.4.3.  

You can install CaSpER R package using the following R commands:

Install dependencies first:

``` r
You may need to install libcurl-devel, libopenssl-devel and libxml2-devel
ex: sudo yum -y install libxml2-devel.x86_64

source("https://bioconductor.org/biocLite.R")
biocLite(c('HMMcopy', 'GenomeGraphs', 'biomaRt', 'limma', 'GO.db', 'org.Hs.eg.db', 'GOstats'))
```

Install CaSpER R package:
``` r
require(devtools)
install_github("akdess/CaSpER")
```

For extracting B-allele freqeuncies from RNA-Seq bam files download BAFExtract C++ source or binary code from [here](https://github.com/akdess/BAFExtract). 


After downloading  source code type the following: 
``` r
cd BAFExtract
make clean
make

```
The executable is located under directory /bin. 


Usage
----------


Step 1. Merge single cell RNA-Seq bam files (required for single-cell studies)

```{bash} 
	 bamtools merge -list <bam_file_names> -out <sample_name>_merged.bam 
	 samtools index <sample_name>_merged.bam
```

Step 2. Extract BAF values from RNA-Seq bam files
	
```{bash}
	samtools view <bam_file> | ./BAFExtract -generate_compressed_pileup_per_SAM stdin <genome_list> <sample_dir> 50 0; ./BAFExtract -get_SNVs_per_pileup  <genome_list> <sample_dir> <genome_fasta_pileup_dir> 20 4 0.1 <output_baf_file>
```
<sample_dir>: the name of sample directory
<output_baf_file>: final output

You can download and unzip generate genome_fasta_pileup_dir files from : 

[for hg38](https://www.dropbox.com/s/ysrcfcnk7z8gyit/hg38.zip?dl=0)

[for hg19](https://www.dropbox.com/s/a3u8f2f8ufm5wdj/hg19.zip?dl=0)
	
To generate genome_list file: 

first download http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes

for hg38:
```{bash}
	fetchChromSizes hg38 > hg38.list
```
	
for hg19: 
```{bash}
	fetchChromSizes hg19 > hg19.list
```

Step 3. Run CaSpER R package (See tutorials below)


Tutorials
----------

1. [Yale meningioma Bulk RNA-Seq dataset](/demo/meningioma.R)
2. [TCGA-GBM Bulk RNA-Seq dataset](/demo/tcga_GBM.R)
3. [GBM Single-cell RNA-Seq dataset](/demo/sCellGBM.R)