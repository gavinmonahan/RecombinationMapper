# Recombination Mapper

Python package for recombination rate estimation and haplotype block and structural variation detection.
**Manuscript in preparation


## Installation

Easiest way to install everything is in a conda environment, required packages are listed in requirements.txt.

1.	`git clone https://github.com/gavinmonahan/RecombinationMapper.git`
2.	`conda create -n recombination_mapper --file ./RecombinationMapper/requirements.txt`


## Usage

A wrapper bash script can be editted which calls all of the python scripts.
  `bash recomb_mapper.sh`
  
For the most part, the scripts will find the relavent file automatically. The only required parameters are the VCF file and height/groups (see `cluster.py`), therefore you should run cluster.py once to find an appropriate height/groups value before running the wrapper script. Alternatively, the individual python scripts can be called separetely.

sci-kit allel is used for parsing VCF files and header errors can cause problems. Escpecially where the Genotype field is not defined in the header.

### cluster.py

Cluster.py clusters individuals together (as the name would suggest) using Ward clustering. Two clustermaps; before and after clustering are exported as well as a dendrogram. Users can specify the height they wish to cut the dendrogram (`-h --height`) or the number of groups they wish to form (`-n --num_groups`). This results in a population groupings file which contains each individuals name and the group they belong to, which can be used for introgression analysis (D-suite or ABBABABAwindows) and for grouping in recomb_rates.py (see below).

```
usage: cluster.py [-h] [-v VCF] [-y HEIGHT] [-n GROUPS] [-o OUTFOLDER]

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     Filtered VCF File (required)
  -y HEIGHT, --height HEIGHT
                        Height to cut dendrogram
  -n GROUPS, --num_groups GROUPS
                        Number of groups desired
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Folder name for output
```

### allele_blocks.py

Allele_blocks.py creates allele blocks (continuous runs of alleles, see manuscript) given a filtered VCF input. The `--phased` flag should only be used if all genotypes are phased.

```
usage: allele_blocks.py [-h] [-v VCF] [-o OUTFOLDER] [-p]

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     Filtered VCF File
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Folder name for output
  -p, --phased          If genotypes are ALL phased
  ```

### merge_blocks.py

As the name suggests, merges adjacent allele blocks where the first block end is within n snps and l length of the start of the next block.

```
usage: merge_blocks.py [-h] [-t TXT] [-o OUTFOLDER] [-n] [-l]

optional arguments:
  -h, --help            show this help message and exit
  -t TXT, --txt TXT     TSV file containing blocks generated by recombmapper.py
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Folder name for output
  -n, --numsnps         Maximum distance for merging two blocks (SNPs). Default = 50
  -l, --length          Maximum distance for merging two blocks (bp). Default = 2000000
  ```

### recomb_rates.py

Reomb_rates.py is the final script and has the function for calculating variable recombination rates. It also calculates SNP density and plots recombination rates and allele blocks grouped using the population groupings file from `cluster.py`. At a minumum, only `allele_blocks.py` and `recomb_rates.py` are needed to calculate recombination rates.

```
usage: recomb_rates.py [-h] [-t TXT] [-t2 TXT2] [-v VCF] [-b {25..100000}]
                       [-s SETS] [-c CENT] [-m MISS] [-n {2..1000}]
                       [-l {1000..10000000}] [-pb] [-pr] [-o OUTFOLDER]

optional arguments:
  -h, --help            show this help message and exit
  -t TXT, --txt TXT     TSV file containing merged blocks generated by
                        recombreduce.py
  -t2 TXT2, --txt2 TXT2
                        TSV file containing blocks generated by
                        recombmapper.py
  -v VCF, --vcf VCF     Filtered VCF file for SNP density
  -b {25..100000}, --bins {25..100000}
                        Number of bins for plotting recombination rates and/or
                        SNP density. Default = 100.
  -s SETS, --sets SETS  Sets file containing individual grouping
  -c CENT, --centromeres CENT
                        TSV file containing centromere locations
  -m MISS, --missing MISS
                        TSV file containing missingness (from allele_blocks.py)
  -n {2..1000}, --numsnps {2..1000}
                        Minimum number of SNPs required to define a block. Default = 2.
  -l {1000..10000000}, --length {1000..10000000}
                        Minimum length (bp) required to define a block. Default = 1000.
  -pb, --plotbar        Plot bargraphs
  -pr, --plotrecomb     Plot recombination rates
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Folder name for output
```
