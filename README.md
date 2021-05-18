# RecombinationMapper

Python package for recombination rate, haplotype block, and structural variation estimation


# Installation

Easiest way to install everything is in a conda environment, required packages are listed in requirements.txt.

1.	git clone https://github.com/gavinmonahan/RecombinationMapper.git
2.	conda create -n recombination_mapper --file ./RecombinationMapper/requirements.txt


# Usage

A wrapper bash script can be editted which calls all of the python scripts.
  `bash run_all.sh`
  
Each script can be called separetly

### cluster.py

Cluster.py clusters individuals together (as the name would suggest) using Ward clustering. Two clustermaps; before and after clustering are exported as well as a dendrogram. Users can specify the height they wish to cut the dendrogram (`-h --height`) or the number of groups they wish to form (`-n --num_groups`). This results in a population groupings file which contains each individuals name and the group they belong to, which can be used for introgression analysis (D-suite or ABBABABAwindows) and for grouping in recomb_rates.py (see below).

```
usage: cluster.py [-h] [-v VCF] [-y HEIGHT] [-n GROUPS] [-o OUTFOLDER]

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     Filtered VCF File
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

### recomb_rates.py

Reomb_rates.py is the final script and has the function for calculating variable recombination rates. It also calculates SNP density and plots recombination rates and allele blocks grouped using the population groupings file from `cluster.py`. At a minumum, only `allele_blocks.py` and `recomb_rates.py` are needed to calculate recombination rates.

