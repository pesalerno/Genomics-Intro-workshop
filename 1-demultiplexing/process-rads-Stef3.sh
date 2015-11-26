#!/bin/bash
#PBS -q q8 1w
#PBS -M patriciasalerno@gmail.com
#PBS -m abe
#PBS -N test-e2



module load stacks/1.35

mkdir ./Pati/raw-data/Lib-1/process_rads_1e2

process_radtags -P -p ./Pati/raw-data/Lib-1/raw/ -b ./Pati/raw-data/Lib-1/barcodes-Stef3-e.txt -o ./Pati/raw-data/Lib-1/process_rads_1e2/  -c -q -r -D --inline_index --renz_1 sphI --renz_2 mspI -i gzfastq 

