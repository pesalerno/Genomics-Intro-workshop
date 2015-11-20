
#Genomics bioinformatics intro workshop


The following document follows a step-by-step protocol for processing genomic data, starting from downloading the dataset into your cluster space and finishing with some simple population analyses in R. It assumes little to no knowledge of unix, servers, and programming.


###1. Getting started
/ / / / / DISCLAIMER: I've never taken a bioinformatics course, it's all been self-taught and picked up in little pieces here and there. There are for sure better ways for doing all of these things, but my aim with this document is to be clear and helpful for people who are starting out just like I did, since I never really had an easy time starting due to being self-taught and not "a natural" at coding. However, the main lesson I've learned so far: ***learn how to google!!!*** everything you need to know is out there, it's a matter of knowing enough of the language to know how to search and what to search./
/
/
/
/
/


Before anything, make a directory for the project, in this case we are calling it Genomics-Bioinformatics-2015. One tip for coding in unix, always make names of files and directories have ZERO spaces in them.... ZERO. 

	mkdir Genomics-Bioinformatics-2015

---> tip: if you have long names that start with unique identifiers, you can always just type the first letter and *"tab complete"* the rest of the file name.

Move into that directory

	cd Genomics-Bioinformatics-2015  ###try the tab completion trick for this command



Second, clone my git repository for the course. Even though you can do this easily on the web and without an account, let's try doing it on terminal! Github is like the GenBank of coding.... let's give it a try, since it's really good practice to share code, to collaborate efficiently through version control software, and to start on this as early as possible!! (I certainly started too late).

1.Create a git account online [here](https://github.com/join), choosing the free plan and open repositories.

2.Download/Install git on your computer [here](http://git-scm.com/downloads)
		

3.Initialize and fork/clone repository for the class into the directory for the course
 
	git init ##within directory that you want the repository to be cloned to
	git clone https://github.com/pesalerno/Genomics-Intro-workshop

Now you should have all the materials for the course organized in the different directories. 

#####Side bar: more on git 
If you want to practice more with git, now that you have an account and the program installed, you can make your new repository (***do this online instead of locally to simplify things***) for your own files/changes that you have for this workshop, within terminal in the following way:

	git init ##within the folder you want to work from 
	git clone <web-address-of-new-directory>
---> make whatever changes, add files, folders, create new text documents, then if you type:

	git status
	
you should see all the changes you've made to your repository. You can add the individual files or all at the same time, with:

	git add <name-of-file.txt>

if you check the status again, you will see this one file changed from "modified" to "added and waiting for commit", then you can commit to the change you just made:

	git commit -m "add informative message for this change"

and finally, once you've made your commit(s) to your change(s), then you can "push" back these changes to the original online repository (so far nothing has been changed other than on your personal computer):

	git push origin master

and now if you refresh the page you should have the new changes already updated on your online repository!!


###2. Downloading files into the raw data directory

First, you need to know where the file is located in the web (i.e. directly from the sequencing facility), then you can copy/paste the address into the terminal window **on the cluster** (from the directory where you wish to copy to) and use the wget command for downloading as follows:

	mkdir raw-data ##make new directory for your data
	cd raw-data
	wget <<web-address>>
	ls ##after downloading is done, do you see all the files you expected to see?

Alternatively, if moving raw data files from your own computer to the cluster, then you simply secure copy the file, which prompts for a password:

	pwd ##on cluster, get the current working directory for the raw data
	scp raw-data-file.gz person@clustername.institution/path/to/file/destination


###2. Adding the path for stacks in your bash profile 

The server needs to know which paths (directories) to look for when you run a program. Thus, you need to set up paths for specific programs that are not directly within your own home directory in the cluster/server. This file is hidden, but you can see it (and whatever is in it) by typing:

	cat .bash_profile

This command simply prints the document, but you cannot alter it. The reason it's invisible is so that it's a bit protected, because your paths are important!! To edit this file, write:

	nano .bash_profile
	
Now you can edit it. Add in a new line to your bash profile as follows:

	PATH=$PATH:/opt/software/stacks-1.26/bin/

Now it will know to search within that directory when you give it a given Stacks command. 


###3. Demultiplexing your dataset with Stacks.

NOTE: **The first step in any RADseq library is always demultiplexing. You are now picking out the barcode reads which are found at the beginning of your illumina read in order to separate them into the individual samples. You only need two things, the code that you will use, and a barcodes/samples file where you have BOTH barcodes (adapter and primer index) and the name you want you sample to be (otherwise the file name will be the barcode, which is zero informative for any human being). Naming the files in a smart way will save headaches down the line. Also, if you have repeated barcodes across different libraries, then not changing the names from the default barcode names will mix your samples eventually. ** 

If you have separate libraries with overlapping barcodes, you need to demultiplex them separately since that is the only sample identifier you have (if you have combinatorial barcodes as in ddRAD, then as long as both adapter and PCR index don't overlap they can be demultiplexed together). 

De-multiplexing is performed using the program [process_radtags](http://creskolab.uoregon.edu/stacks/comp/process_radtags.php) individually for each library within its directory, and renamed with sample names within Stacks using the appropriate barcodes/names text files, found [here](https://github.com/pesalerno/Genomics-Intro-workshop/tree/master/1-demultiplexing). The barcodes file is a simple text-delimited file with first column being adapter, second column being primer index, and third being the final file name you want (ideally an individual sample code, and maybe a locality code as well). 

**WARNING**: _NEVER_ edit text files in word or something similar... it needs to be in a simple text editor such as text wrangler or BBedit. Also, always edit these files while seeing "invisible characters". 

You need to have the appropriate barcode files within the appropriate library folder if demultiplexing libraries separately. 

Figure out how your barcodes are set up within your sequence file, in order to determine how to set up the process_radtags code.

To look into a gzipped file (for example, with head):

	zcat file.name.fastq.gz | head ##may need to do gzcat instead

Or you can also do less:

	zless ##or gzless

	

The commands for process_radtags for the first single-end library that we will analyze are:

	process_radtags  -b barcodes-a -o ./process_rads_B/  -q -D -w 0.15 -s 10 
		--inline_index --renz_1 sphI --renz_2 mspI -f ./Stef_3_ATCACG_L008_R1_001.fastq.gz -i gzfastq 

	(example for library #1610 for *Xantusia*)

Process radtags also cleans your data as it demultiplexes. 
 
This command needs to be set up within the job scheduler (file of type ***.sh***), so set it up in the cluster (either set up the file on your computer and scp to cluster, or edit in cluster with ***nano*** text editor).

To start a text file from scratch on the cluster (in general in unix) type:

	touch process-rads-A.sh

then open the blank file to beign editing:


	nano process-rads-A.sh

copy/paste the text for the job scheduler on the first line of the file, then copy/paste the commands for process_radtags. The file should be in the directory where your raw data folder is located (not within it) and then you can specify where your raw reads are located with the raw reads folder name. 

To make the other process_radtags shell script for submitting the second demultiplexing job, copy the file into the same directory in order to save with a different name:

	cp process-rads-A.sh process-rads-B.sh

now edit the second .sh file to make sure it's running the appropriate data and the bleh bleh bleh. 


 
#####3b. Merge libraries into single ***denovo*** directory

Copy all renamed libraries for all individuals into the ***denovo*** directory; after checking that copying was performed successfully, then delete the files in the previous directory. 

	cp *.fq /path/to/file   ##this will copy all files that end in .fq into the given directory 
								(asterisk is a the wildcard)
	cd /path/to/file
	ls -ltr     ##do files look ok and of the right size? 
	rm -r directoryname    ##remove the entire directory were the demultiplexed sequences were originally 
								demultiplexed (or you can also just delete the files that ended 
								up being copied over only)
	

#####3c. Looking at our demultiplexed dataset

Now that we've separated our dataset into our individuals, we can see how much data we have per individual. The fastest but not ideal way to look at this is simply by looking at the sizes of the .fq files that we now have for each individual. To get a list with file size information, type:

	ls -l

The files are huge!! and it's hard to read. So how can we look at them in a different unit (Megabyte? Terabyte?). We can find out by looking at the manual:

	man ls

(tip, scroll down to -h). The you'll se that we can use that flag:

	ls -l -h

and now we see mor easily the sizes of our files. But can we order them? Let's pipe our list command with sort, using -rn (reverse numeric) to order from larger to smaller file:

	ls -l -h | sort -k 5 -rn

However, if we want a more precise piece of information, we need read counts instead! 

	echo -e 'SAMPLE_ID_FULL\tNUM_READS'
		for file in *.fq
	do
    	echo -n $(basename $file .fq)$'\t'
    	cat $file | grep '^@.*' | wc -l
	done

Here the syntax is cmore complicated, but in essence we are "grepping" (finding and grabbing) the first characters of each line of sequence, so in essence simply counting the lines of sequence within each file. 
	
###4. Running denovo map permutations for parameter testing


Running [denovo_map](http://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php) is not a complicated process... but since different datasets can be very sensitive to parameter settings, we need to explore our dataset with different parameter combinations in order to decide what's the best approach (keeping good data, but not overfiltering). 

Denovo-map needs to be run with all your data, and even though it can be threaded for faster computation speed, the more individuals you have, the longer it takes, since the last step can't be threaded (all stacks need to be compared at once to build the catalog loci).

The main parameters to vary are:

m — specify a minimum number of identical, raw reads required to create a stack.

M — specify the number of mismatches allowed between loci when processing a single individual (default 2).

n — specify the number of mismatches allowed between loci when building the catalog (default 0).


**Note 1**: The higher the coverage, the higher the m parameter can be. 

**Note 2**: M should not be 1 (diploid data) but also should not be very high since it will begin to stack paralogs. 

**Note 3**: n will depend on how divergent our individuals/populations are. IT should not be zero, since that would essentially create zero SNPs, but 1 also seems unrealistically low (only a single difference between individuals in any given locus), so in these kinds of datasets we should start permutations starting from 2.  If you use n 1 it is likely to oversplit loci among populations that are more divergent. 



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


In order to make the input file for denovo_map, you need to generate a popmap, with all the correct file names form the stacks output, thus it's easier to just save a list of the .fq files into a text document to start that process. 

	ls *.fq > sequences-list.txt

The download that document on to your computer (from the shell window that is NOT on the cluster but on your own machine):

	scp user@cluster.name:/path/to/file/sequences-list.txt . 

the dot at the end of that command tells the address of where you want the file saved (in this case, it's just the current directory you're in). 

You now open that file in a text editor, and through find/replace commands (essentially using **regular expressions** but with your text editor) you transform your list into the input file for denovo_map, which needs EVERY SINGLE input file listed on the .sh file and it's path. Each line should look like this:

	-s ./IND-a-sequence.fq





----------------------------------------------
----------------------------------------------




/

/

/

/

/

/

/

/

/






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

 