
#Genomics bioinformatics intro workshop


The following document follows a step-by-step protocol for processing genomic data, starting from downloading the dataset into your cluster space and finishing with some simple population analyses in R. It assumes little to no knowledge of unix, servers, and programming.


###1. Downloading files into the raw data directory
Before anything, make a directory for the project, in this case we are calling it Genomics-Bioinformatics-2015, so:

	mkdir Genomics-Bioinformatics-2015

****tip: if you have long names but that start with unique identifiers, you can always just type the first letter and *"tab complete"* the rest of the file name.

Move into that directory

	cd Genomics-Bioinformatics-2015  ###try the tab completion trick for this command


Second, you need to know where the file is located in the web (i.e. directly from the sequencing facility), then you can copy/paste the address into the terminal window (from the directory where you wish to copy to) and use the wget command for downloading as follows:

	mkdir raw data ##make new directory for your data
	cd raw data
	wget <<web-address>>
	ls ##after downloading is done, do you see all the files you expected to see?

Alternatively, if moving raw data files from your own computer to the cluster, then you simply secure copy the file, which prompts for a password:

	pwd ##on cluster, get the current working directory for the raw data
	scp raw-data-file.gz person@clustername.institution/path/to/file/destination


###2. Making sure your path is set up properly for the program to run and the data to be found

I'm not sure what to write here, need to look this up a bit better, but certainly something we should check before getting started. 


###3. Demultiplexing your dataset.

If you have separate libraries with overlapping barcodes, you need to demultiplex them separately since that is the only sample identifier you have (if you have combinatorial barcodes as in ddRAD, then as long as both adapter and PCR index don't overlap they can be demultiplexed together). 

De-multiplexing will be performed using the program [process_radtags](http://creskolab.uoregon.edu/stacks/comp/process_radtags.php) individually for each library within its directory, and renamed with sample names within Stacks using the appropriate barcodes/names text files, found [here](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/barcodes-1933.txt).

- The commands for process_radtags for the Paired-end libraries are:

		process_radtags -P -p ./PE-lib-1610/ -b barcodes-1610.txt -i gzfastq - \
		o ./processed-1610/ -e sbfI -c -q -r -D
	(example for library #1610 for *Xantusia*)
 


#####3b. Merge libraries into single ***denovo*** directory

Copy all renamed libraries for all individuals into the ***denovo*** directory; after checking that copying was performed successfully, then delete the files in the previous directory. 

	cp *.fq /path/to/file   ##this will copy all files that end in .fq into the given directory 
								(asterisk is called the wildcard)
	cd /path/to/file
	ls -ltr     ##do files look ok and of the right size? 
	rm -r directoryname    ##remove the entire directory were sequences were originally 
								demultiplexed (or you can also just delete the files that ended 
								up being copied over only)
	
	
###4. Running denovo map permutations for parameter testing


It's important to set up permutations to explore the data and how well the parameters are recovering "good" SNPs that have decent coverage within and among populations. Some parameters such as xx and xx are very sensitive to y and z, so bleh. 

####keep editing from here down
Because we have decently high coverage, I'm varying the -m parameter (# reads required to form a stack) from 3-7, skipping even numbers (just to have a bigger range initially). Also, I'm not at all using -n 0 or 1, since they are biologically unrealistic given these datasets, and likely to eliminate and oversplit loci among individuals/populations with higher divergence. 

Permutations | -m | -M | -n | --max_locus_stacks 
------------ | ------------- | ------------ | ------------- | ------------ |
a | 3 | 2 | 2 | 3 | 
b | 5 | 2 | 2 | 3 |
c | 7 | 2 | 2 | 3 | 
d | 3 | 3 | 2 | 3 |
e | 3 | 4 | 2 | 3 |
f | 3 | 5 | 2 | 3 |
g | 3 | 2 | 3 | 3 |
h | 3 | 2 | 4 | 3 |
i | 3 | 2 | 5 | 3 |
j | 3 | 2 | 2 | 4 |
k | 3 | 2 | 2 | 5 |

----------------------------------------------
----------------------------------------------


###Intermediate step: get a quick SNP matrix only with SR sequences for Chris' report

######1. check to see which files are repeats in different libraries, or else when moving them to do denovo_map they will be rewritten/lost.

Libraries #1612 and #1835 are essentially duplicates, with few exceptions. Library #1834 is mostly unique with a few that are repeats from #1612. So for ***Xantusia***, I'm renaming all sequence files and adding a -1.fq.gz to library #1 (1612), -2.fq.gz to library 2 (1834), -3.fq.gz to library 3 (1835), -4.fq.gz to library 4 (1994), and -5.fq.gz to library 5 (1995). Renaming all files at once using the following code:

	rename .fq.gz -1.fq.gz *.fq.gz
	
	


------------------------------------


######2. Merge fasta files for library duplicates

After being renamed, move all files back to the SR-denovo-prelim folder and there I merge the fasta files. I merge following these guidelines (from this [source](http://www.researchgate.net/post/How_do_I_merge_several_multisequence-fasta_files_to_create_one_tree_for_subsequent_Unifrac_analysis)):

*To merge several files use the SHELL, go to your folder where the files are and use the cat command. E.g. to merge seqfile001.fasta, seqfile002.fasta and seqfile003.fasta type*

	cat seqfile001.fasta seqfile002.fasta seqfile003.fasta > seqcombined.fasta


*or if you have more files use*

	cat *.fasta > seqcombined.fasta


The duplicated files are sorted into a separate folder before merging, just to keep track of what's being merged. Then the post-merged files are sorted back into the general directory containing all sequences. Total number of files before merging duplicates from different ***Xantusia*** library preps was 187, and after merging duplicate individuals we now have 142 files for denovo_map input. Total number of files before merging duplicates from different ***Pseudacris*** library preps was 180, and after merging duplicate individuals we now have 132 files for denovo_map input. 

------------------------------------


######3. Running denovo_map for SR reads

The code used for running denovo_map for only the SR reads was:

	denovo_map.pl -m 3 -M 2 -n 1 -T 16 -b 1 -t -S -o ./denovo-1/ \



------------------------------------


######4. running program populations for exporting SNP matrix

I'm running populations with the filters for keeping 50% of SNPs that are present in each populations and with two settings for numbers of populations kept. Here is the popmap for [*Xantusia*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/popmap_Xari.txt) and the popmap for [*Pseudacris*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/popmap-Pseu.txt). 

The script for populations for ***Xantusia*** and for **Pseudacris** was:

	populations -b 1 -P ./denovo-1 -M ./popmap_Xari.txt  -t 36 -p 6 -r 0.5 --write_random_snp --structure --genepop --vcf

I ran this twice, first with 6/7 populations (-p 6), and second with all populations (-p 7).I tried to run the filter of minor allele frequency (**--min_maf**) and it kept failing, maybe I need to upgrade version of Stacks, or maybe I don't know how to set it! This is why I saved it in **--vcf** format so that it can be viewed and filtered in [GATK](https://www.broadinstitute.org/gatk/) easily (or other variant calling format software).

The total number of SNPs were:


- ***Xantusia*** p=6 => 7739 snps
	
- ***Xantusia*** p=7 => 1319 snps 

- ***Pseudacris*** p=6 => 36921 snps

- ***Pseudacris*** p=7 => 36921 snps (did I make a mistake in the previous one? this one's correct according to script). it's just odd they're identical.... 


Just out of curiosity, running ***Pseudacris*** again with a much more stringent setting of r=0.8 to see how many SNPs are still kept. 

Based on number of SNPs, it's already suggesting that ***Xantusia*** has a much higher divergence among populations.... 

-> When doing -p7 -r 0.7 for ***Xantusia*** you get ZERO SNPs.... which means that there is either a huge amount of missing data, or that the mainland population is too divergent for this constraint.... 




----------------------------------------------
----------------------------------------------
----------------------------------------------
----------------------------------------------



#

#

##Step 2: create reference genome with PE reads

######2.1. Run denovo_map program for paired-end reads

Run in Stacks the script [denovo_map](http://creskolab.uoregon.edu/stacks/comp/denovo_map.php) with libraries with higher depth of coverage and paired-end reads for creating the reference genome. I'm trying to be highly conservative about the flag parameters. 

	> denovo_map.pl -m 3 -M 1 -n 0 -T 16 -b 1 -t -S -o ./denovo-1/ \


From Hohenlohe et al. 2013: *"We grouped the forward and reverse reads from all individuals in these populations into a separate file for each RAD locus, using the STACKS program sort_read_pairs.pl"*

More from Hohenlohe et al. 2013: *"We created a catalog of RAD tag loci using cstacks and matched indi- viduals against the catalog using sstacks. We populated and indexed a MYSQL database of loci using load_ radtags.pl and index_radtags.pl and then exported the data using export_sql.pl. Finally, we grouped the for- ward and reverse reads from each individual corre- sponding to each RAD locus using sort_read_pairs.pl."*

More from Hohenlohe et al. 2013: *"We assembled the reads in each file separately to produce a set of RAD contigs (Fig. 2b), using both VELVET (Zerbino & Bir- ney 2008) and CAP3 (Huang & Madan 1999) assembly software."*


#

#

#

#

#

#

#

#

#




---> May want to run [exec_velvet.pl](http://catchenlab.life.illinois.edu/stacks/comp/exec_velvet.php) to generate collated fasta file for reference genome.

Use the output consensus sequence file from denovo_map (*"catalogs.tags.tsv"*) and transform to fasta format for input into [bwa](http://bio-bwa.sourceforge.net/bwa.shtml), using Kelly's R script (need to modify once I run it):

	
	# Import the data and check the structure
	tags<-read.table('batch_1.catalog.tags.tsv', header=FALSE)
	# all the sequences are in $V9

	unique(tags$V2) # sample ID is meaningless here
	length(unique(tags$V3)) #417153 -- each row has a unique locus ID

	# verify that all of the tags are consensus:
	consensus.tags<-subset(tags, V6=='consensus')
	length(consensus.tags[,1]) #417153
	length(tags[,1]) #417153

	# each sequence needs a fasta header to uniquely identify it
	fa.id<-paste('>', tags$V3, '_pseudoreference_pe_concatenated_without_rev_complement', sep='')

	fa<-cbind(fa.id, as.character(tags$V9))
	write.table(fa, file='D_variabilis_denovo_psuedoreference.fa', quote=FALSE, sep='\n', row.names=FALSE,
            col.names=FALSE)
 

######2.2. Index reference genome and align in bwa

Create reference genome from fasta consensus sequences:

	> bwa index paired-end.fasta

Then align all of the de-multiplexed files (single-end 100bp reads) to reference genome:

	> bwa mem -t 6 paired-end.fasta all-reads-demultiplexed.fq > align-allRADs.sam
	
Transform .sam file to .bam file for visualization in IGV:

	> module load samtools
	> samtools view -b -S -o 454-align.bam 454-align.sam
	> samtools sort 454-align.bam 454-align.sorted
	> samtools index 454-align.sorted.bam

Visualize in IGV

Use either .sam or .bam alignment for input in ref_map.pl

##Step 3: Map to new reference genome in Stacks

Use ref_map.pl for mapping the raw read files to the reference genome we created using the in-depth paired-end reads plus the single-end reads of all individuals.

	> ref_map.pl -o /path/for/output -n 2 -m 2 -T 12 -O popmap.txt -b 1 -S -s ./all-sequences-here\

##Step 4: Filter data with haplotype corrections in Stacks

Here we essentially re-run the stacks pipeline (post-demultiplexing) to make corrections based on population haplotype statistics (reduced probability of obtainind fake haplotypes created by paralog stacks and PCR duplicate errors)

	> rxstacks -b 1 -P ./input-stacksoutput/ -o /new-output-path/ --prune_haplo --lnl_filter --lnl_dist 

 