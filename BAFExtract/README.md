# BAFExtract

BAFExtract generates the B-allele frequency shifts from mapped reads alone without the need for a heterozygous variant call set. BAFExtract a part of CaSpER (https://github.com/akdess/CaSpER) 

There are two steps:
1. Generation of the pileup.
2. Generation of the B-allele frequency profile.

Please open a new "issue" on github or contact authors Akdes.Harmanci@uth.tmc.edu or arif.o.harmanci@uth.tmc.edu for questions.

# Installation

Type make to build BAFExtract. The executable is built under bin directory. The code is tested on various Unix based systems.

# Usage 

Extract BAF values from RNA-Seq bam files
	
```{bash}
	samtools view <bam_file> | ./BAFExtract -generate_compressed_pileup_per_SAM stdin <genome_list> <sample_dir> 50 0; ./BAFExtract -get_SNVs_per_pileup  <genome_list> <sample_dir> <genome_fasta_pileup_dir> 20 4 0.1 <output_baf_file>
```
<sample_dir>: the name of sample directory
<output_baf_file>: final output

You can download and unzip genome_fasta_pileup_dir files from : 

[for hg38](https://www.dropbox.com/s/ysrcfcnk7z8gyit/hg38.zip?dl=0)

[for hg19](https://www.dropbox.com/s/a3u8f2f8ufm5wdj/hg19.zip?dl=0)
	
You can download genome_list files from : 

[for hg38](https://www.dropbox.com/s/rq7v67tiou1qwwg/hg38.list?dl=0)
	
[for hg19](https://www.dropbox.com/s/jcmt23nmuzm6poz/hg19.list?dl=0) 

# Example
[download example bam file](https://www.dropbox.com/s/1vl6iip0b8jwu66/SRR1295366.sorted.bam?dl=0)

```{bash} 
mkdir test; samtools view SRR1295366.sorted.bam | ./bin/BAFExtract -generate_compressed_pileup_per_SAM stdin hg38.list test 50 0; ./bin/BAFExtract -get_SNVs_per_pileup hg38.list test ./hg38/ 20 4 0.1 test.baf
```
