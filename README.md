# Low_pass
Finding CNAs in low pass WGS data. Use `Rscript low_pass_cna_tool.R --help` to get started. Requires some package installation to get started. 

A docker container with packages loaded can be found here: https://hub.docker.com/r/taytayf/low_pass_cna_bin_replacer

This R tool builds on the bioconductor package QDNAseq in order to generate multiple analyses of whole genome sequencing data. This script runs QDNAseq two times, once for bins of a large size and once for a smaller increment. The script then takes the two segments files and two binned copy number files, and replaces the larger bins with the smaller bins along a region of interest (typically, a focal copy number variant, which only occurs in a region of a few megabases instead of an entire chromosome arm, as in a normal CNV).

## Install
This script requires the following packages:
From bioconductor:
- Biobase
- BSgenome
- genome.Hsapiens.UCSC.hg19
- QDNAseq
- QDNAseq.hg19
From CRAN:
- dplyr
- future
- ggplot2

Instructions for installing and using the docker container are on the docker hub, linked above.

## Usage
The script can be run using `Rscript low_pass_cna_tool.R`. `Rscript --vanilla` can be specified to ignore R environments.  

Use `Rscript low_pass_cna_tool.R --help` for more details.  

### Mandatory arguments
`-b` or `--bamfile=FILE` specifies the .BAM file with which to run the tool on.  

Two options exist for specifying a region of interest, `-g --gene=GENE` or `--chromosome`, `--start`, and `--end`. Specifying a gene will create a 150 kilobase pair region surrounding that gene in which to look for amplifications. If `-g` is not used, chromosome, start, and end _must_ be specified.  

### Optional arguments
`--smallbins` and `--largebins` control the size of the bins to be used. There are a finite number of bins that can be used: 1, 5, 10, 15, 30, 50, 100, 500, 1000. The large bin size must be a multiple of the small bin sizes. If not specified, small bin size will be 50 kilobase pairs and large bin size will be 500 kilobase pairs.

`--threads` controls the number of compute threads that can be used for the QDNAseq steps.  

`--zerosegments` changes the default QDNAseq segmenting function to include segments where copy number is approximately 0, which is excluded by default.  

`--saveplots` saves some plots that are output by QDNAseq, as well as a plot describing the number of bins in the region of interest.  

`--warnings` enables or disables warnings for number of bins and reads within them in certain area.  

`--out` specifies a character string to be used as a part of the output file names.

## Further Development
For the time being, only a handful of genes are specified by name, as the data to determine which segments were applicable was limitted. If someone, in the future, wished to change which genes were available and their associated regions, they would just need to change the `data.frame` that begins on line 11 in `low_pass_cna_tool.R`.  

The tool also uses a t-test to determine statistical differences between bin regions, which could be adjusted under the section `#statistics`.

There is another script, `cpdm_sample_processing.R`, which uses a few input `.seg` and `.log2ratio` files in order to define the regions used in each gene, as well as the amplification level threshold.  

