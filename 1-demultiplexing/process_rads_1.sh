#!/bin/bash 
#
#SBATCH -J process_rads_TEP
#SBATCH -o process_rads_TEP.%j.out
#SBATCH -e process_rads_TEP.%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#
#SBATCH --mail-user=patriciasalerno@gmail.com
#SBATCH --mail-type=ALL
#
#SBATCH -A Frog-RadSeq



module load stacks


process_radtags -p ./raw_data_A/ -b adapters_Tep -o ./process_rads_A/  -c -q -r -D --inline_index --renz_1 sphI --renz_2 mspI -i gzfastq 
	

