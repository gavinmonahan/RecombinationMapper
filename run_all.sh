
#!/bin/bash

echo -e "\nWelcome to Recombination Mapper\n"

output="recomboutput"
VCF="tamerged01.vcf.recode.vcf"
centromeres="../../triticumaestivumcentromeres.tsv"

python3 cluster.py -v ${VCF} -y 310 -o ${output}
python3 allele_blocks.py -v ${VCF} -o ${output}
python3 merge_blocks.py -o ${output}
python3 recomb_rates.py -v ${VCF} -c ${centromeres} -o ${output}
