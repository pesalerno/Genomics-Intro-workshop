module load stacks/1.35

mkdir ./Pati/raw-data/Lib-1/process_rads_1e2

process_radtags -P -p ./path/to/raw/data/ -b ./path/to/folder/barcodes-Stef3-e.txt -o ./path/to/folder/process_rads_1e2/  -c -q -r -D --inline_index --renz_1 sphI --renz_2 mspI -i gzfastq 

