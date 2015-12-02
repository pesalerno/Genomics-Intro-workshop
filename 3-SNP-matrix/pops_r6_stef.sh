#!/bin/bash
#PBS -q q8 1w
#PBS -M patriciasalerno@gmail.com
#PBS -m abe
#PBS -N pops-r6

cd $PBS_O_WORKDIR

module load stacks/1.35

populations -b 2 -P ./populations_r_6/ -M ./popmap_Stefania.txt -k -p 3 -r 0.6 -a 0.02 -f p_value -t 36 --structure --genepop --vcf --write_random_snp
