
#Genomics bioinformatics intro workshop


The following document follows a step-by-step protocol for processing genomic data, starting from downloading/uploading the raw data into your cluster space, genotyping with Stacks, and finishing with some simple population analyses in R. It assumes little to no knowledge of unix, servers, and programming.


/ / / / / DISCLAIMER: I've never taken a bioinformatics course, it's all been self-taught and picked up in little pieces here and there. There are for sure better ways for doing all of these things, but my aim with this document is to be clear and helpful for people who are starting out just like I did, since I never really had an easy time starting due to being self-taught and not "a natural" at coding. However, the main lesson I've learned so far: ***learn how to google!!!*** everything you need to know is out there, it's a matter of knowing enough of the language to know how to search and what to search.


Also, [this book](http://www.amazon.com/Practical-Computing-Biologists-Steven-Haddock/dp/0878933913) is one of the best resources you will find out there for starting out and learning the different languages that are most relevant for biologists. 
/
/
/



###1. Getting started

Before anything, make a directory for the project, in this case we are calling it Genomics-Bioinformatics-2015. One tip for coding in unix, always make names of files and directories have ZERO spaces in them.... ZERO. 

	mkdir Genomics-Bioinformatics-2015

---> tip: if you have long names that start with unique identifiers, you can always just type the first letter and *"tab complete"* the rest of the file name.

Move into that directory

	cd Genomics-Bioinformatics-2015  ###try the tab completion trick for this command


#####1.2. Git-ing started (yep, I'm hilarious)

Git is an excellent free online resource that acts like a repository for yourself and for sharing with collaborators. You can pay and make private repositories, but the whole idea is open source coding and sharing. Github is like the GenBank of coding.... let's give it a try, since it's really good practice to share code, to collaborate efficiently through version control software, and to start on this as early as possible!! (I certainly started too late).

Start by cloning my git repository for the course, so that you'll have all the material you need. Even though you can do this easily on the web and without an account, let's try doing it on terminal! 

1.Create a git account online [here](https://github.com/join), choosing the free plan and open repositories.

2.Download/Install git on your computer [here](http://git-scm.com/downloads)
		

3.Initialize and fork/clone repository for the class into the directory for the course
 
	git init ##within directory that you want the repository to be cloned to
	git clone https://github.com/pesalerno/Genomics-Intro-workshop

Now you should have all the materials for the course organized in the different directories. 

#####Side bar: more on git 
If you want to practice more with git, now that you have an account and the program installed, you can make your own new repository (***do this online instead of locally to simplify things***) for your own files/changes that you have for this workshop, within shell in the following way:

	git init ##within the folder you want to work from 
	git clone <web-address-of-new-directory>
---> make whatever changes, add files, folders, create new text documents. For example, you could create a new blank document with:

	touch blank-example.txt

then if you type:

	git status
	
you should see all the changes you've made to your repository. You can add the individual files or all at the same time, with:

	git add <name-of-file.txt>

if you check the status again, you will see this one file changed from "modified" to "added and waiting for commit", then you can commit to the change you just made:

	git commit -m "add informative message for this change"

and finally, once you've made your commit(s) to your change(s), then you can "push" back these changes to the original online repository (so far nothing has been changed other than on your personal computer):

	git push origin master

and now if you refresh the page you should have the new changes already updated on your online repository!!


Note: In order to make new directories in git, the directory will not be added unless it has a document in it (in essence, what you are adding is a new file within that directory, and the directory gets cloned into the repository), so when you are ready to add a file within a new directory, then add/commit/push the change and you should be set. 


#####1.3. Downloading files into the raw data directory

First, you need to know where your data files are located in the web (i.e. directly from the sequencing facility), then you can copy/paste the address into the terminal window **on the cluster** (from the directory where you wish to copy to) and use the wget command for downloading as follows:

	mkdir raw-data ##make new directory for your data
	cd raw-data
	wget <<web-address>>


Alternatively, if moving raw data files from your own computer to the cluster, then you simply secure copy the file, which prompts for a password:

	pwd ##on cluster, get the current working directory for the raw data
	scp raw-data-file.gz person@clustername.institution:/path/to/file/destination

For this specific dataset, we have two "libraries" (four sequence files, one with each barcode), and they are located in my directory under raw-data:

	/lustre/home/ciencias/biologia/pregrado/s.herrera706/Pati/raw-data

So, from your own home directory, make a raw-data file and copy all the raw data files from my directory into yours:

	mkdir raw-data
	cp /lustre/home/ciencias/biologia/pregrado/s.herrera706/Pati/raw-data/Stef* .
	
The DOT at the end is crucial... you're telling the shell **where** you want that file copied to, which in this case is your current directory. Notice that I added an asterisk, which is the wildcard. This will copy ALL files within that directory that start with ***Stef***. Alternatively, you could copy all files of fastq.gz format as such:

	cp /lustre/home/ciencias/biologia/pregrado/s.herrera706/Pati/raw-data/*fastq.gz .



#####1.4. Let's start messing with the raw data.

The example dataset we are using is made of only fifty individuals, multiplexed in a single Illumina lane, obtaining 150bp Paired-end reads. In reality, I have a single lane of another set of 16 individuals, which I've used to generate a reference genome. But, for now, we will keep it simple and only use those two pools of data. 


The names of the raw data files are:

	Stef_4b_CGATGT_L008_R1_001.fastq.gz
	Stef_4b_CGATGT_L008_R2_001.fastq.gz
	Stef_3_ATCACG_L008_R2_001.fastq.gz
	Stef_3_ATCACG_L008_R1_001.fastq.gz

Notice we have two files for each pool, which correspond to each read of our paired-end Illumina data (R1 and R2). The files are formatted as fastq, which is the standard for most sequence files, but the termination is gz, which is a type of file compression (g-zipped). The best way to transfer files and store them is zipped, but sometimes we want to uncompress them, for example to look into the files for some reason... in this case, to uncompress you would simply type:

	gunzip name-of-file.gz

And it expands file to original format. Let's do that to look into one of the files.... but it may take a bit...

Now let's look into the file... normally, for a text file we would do something like

	cat name-of-file.fastq

which prints to screen the entire file in a non-editable format, which is awesome for looking at many text files... HOWEVER.... if you did this with this file.... well, it has millions of lines corresponding to millions of reads!!! (if you did, you can now "quit" the cat command by typing Control C).

A better way of loking into these gigantic files is, for example, head or tail, to look at the beginning and end of files, respectively:

	head name-of-file.fastq 
	
The default is only ten lines, so if we want to see more, we can type:

	head name-of-file.fastq -n 25

Now we can have a look at the file format. Each new sequence begins with @, but the first line of the sequence is some standard Illumina information, which you can find [here](http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm).

If you do the same but with the command tail, you will see the exact same format until the end of file... just lines and lines of sequence. Since we now know that each sequence begins with @, if we want to find out how many total sequences we have, we can simply count those characters.... right? Let's try it:

	grep -c '@D3' Stef_3_ATCACG_L008_R1_001.fastq

This will find and count the number of lines that start with the argument '@D3'. It will take a few seconds.... How many do we have?

Awesome!! now we can get started with the real stuff.... analyzing/processing our data with STACKS! 

**NOTE**: You can also look into a gzipped file with slightly modified commands/programs, for example:

	zless ##or gzless




###2. Demultiplexing your dataset with Stacks.

NOTE: **The first step in any RADseq library is always demultiplexing. You are now picking out the barcode reads which are found at the beginning of your illumina read in order to separate them into the individual samples. You only need two things, the code that you will use, and a barcodes/samples file where you have BOTH barcodes (adapter and primer index) and the name you want your sample to have (otherwise the file name will be the barcode, which is zero informative for any human being). Naming the files in a smart way will save headaches down the line. Also, if you have repeated barcodes across different libraries, then not changing the names from the default barcode names will mix your samples eventually. ** 


#####2.1. Setting up the barcodes files

De-multiplexing is performed using the program [process_radtags](http://creskolab.uoregon.edu/stacks/comp/process_radtags.php) individually for each library within its directory (you can do it all at once, if barcodes are not the same, but let's do it separately for repetition and practice).

Other than the raw data, we need only a single input file for this step, the barcodes file. This file has all the info to pick out the combinatorial barcodes from your raw sequence reads, and it will split them up by sample name, according to your file. I have written up a [single](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/1-demultiplexing/barcodes-Stef-3.txt) barcodes file, so you need to build the other. 

You can secure copy this file that is now on your git directory on your computer into your directory in the cluster:

	scp barcodes-Stef-3.txt username@clustername:path/to/directory



The barcodes file is a simple text-delimited file with first column being the unique adapter, second column being the primer index, and third being the final file name you want (ideally an individual sample code, and maybe a locality code as well, whatever is informative for you later down the pipeline). In this case, we have ***Locality_sampleID*** format for names (we only have three localities for this project). 

Now, build the same input file but for library ***Stef_4***. The names of the sequences are these:

	Er_413
	Er_414
	Er_416
	Er_418
	Er_419
	Er_420
	Er_422
	Er_423
	Er_424
	Er_425
	Er_426
	Er_427
	Er_428
	Er_429
	Er_431
	Er_432
	Er_433
	Er_434
	Er_435
	Er_436
	Er_467
	Er_468
	Er_469
	Er_470
	TNHC05833


And it has the exact same adapters as the first library, but with a different Index primer (can you figure out which one it is??). 

**WARNING**: _NEVER_ edit text files in word/excel or something similar... it needs to be in a simple text editor such as text wrangler or BBedit. Also, always edit these files while seeing "invisible characters"... it's good practice to avoid unnecessary headaches! 

#####2.3. De-multiplexing and commands

Let's demultiplex these two libraries separately for now, to learn a bit more by repeating and adding a couple of steps. To do this, create a de-multiplex directory for each one of the libraries, and also a raw directory to put the raw data in (in each of the library directories) so that steps are carried out cleanly and separately, for now...

You need to have the appropriate barcode file within the appropriate library folders if demultiplexing libraries separately. 

Figure out how your barcodes are set up within your sequence file, in order to determine how to set up the process_radtags code (doing any of the commands we did earlier to look into the files).


Now that you know how your sequence file looks like (it will vary from one facility to another), then you can enter the option for whether the barcode occurs in line with the sequence or not, and the other options as well (see the [process_radtags](http://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) manual for more info).
	

The commands for process_radtags for the first paired-end library that we will analyze can be found [here](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/1-demultiplexing/process-rads-Stef3.sh), but you need to edit the path to your working directory and to your raw data. 

The commands for the job scheduler specific for the Uniandes Cluster can be found [here](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/1-demultiplexing/cluster-header.txt.sh). You need to modify the cluster header with your personal information, and the name you want to give your job. Try to mofidy this file in the shell using the program nano, as such:

	nano cluster-header.txt

Now you can edit the file, and follow the instructions on the bottom of the shell window to save and exit the file. 


We can set up our job file in one of two ways... you can either copy/paste the job scheduler header (specific to each cluster) and after that the Stacks commands below (using the program nano on the cluster), or we can merge the two files back to back using cat. To merge two files back to back using cat, you only need to do:

	cat cluster-header.txt process-rads-Stef3.sh > process-rads-Stef3-b.sh 


To start the text file from scratch on the cluster (in general in unix) so that you can just copy/paste the job scheduler and the Stacks commands, type:

	touch process-rads-Stef3.sh

then open the blank file to beign editing:


	nano process-rads-Stef3.sh

copy/paste the text for the job scheduler on the first line of the file, then copy/paste the commands for process_radtags. The file should be in the directory where your raw data folder is located (not within it) and then you can specify where your raw reads are located with the raw reads folder name. 


To make the other process_radtags shell script for submitting the second demultiplexing job, copy the file into the other demultiplexing directory in order to save with a different name:

	cp process-rads-Stef3.sh ./Stef4/process-rads-Stef4.sh

now edit the second .sh file to make sure it's running the appropriate data and within th appropriate directories. 

#####2.4. Submit job on cluster

To submit your job on the cluster (this may be different depending on cluster), all you have to do is write the submit command followed by the name of your file:

	qsub name-of-job.sh

To see if your job is running properly, type:

	qstat user-name

And you should see the information on the jobname, where the job is running, and the status of the job. To kill a job, you can type:

	qdel jobName

And the jobName is the one you see when you qstat the job you just submitted. 


#####2.5. Look at your logfiles once your job is done
If you set up your cluster header (job scheduler) script correctly, you should recieve an email when you submit the job, and another one when the job is done. When the job is done, you should cd into your job folder and do:

	ls -l 
You will see all the files that have been generated as output. We have the discarded files per each file, and the forward and reverse reads for all of them as well. Now look into your logfile, which should look almost exactly (if not identical) to [this](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/1-demultiplexing/process_radtags.log).

 
#####3b. Merge de-multiplexed libraries into single ***denovo*** directory

Copy all renamed libraries for all individuals into a ***denovo-map*** directory; after checking that copying was performed successfully, then delete the files in the previous directory. 

First clean your directory to facilitate copying the desired files to your denovo-map directory, from within your process-rads directory:

	mkdir discards
	mv *.rem.* discards ###this will move all files that contain '.rem.' into the discards directory (asterisk is the wildcard)
	mv *gz.discards discards

Now your directory only contains the logfile and the "kept" sequences for each individual. Now move your files to denovo-map:

	cp *.fq.gz /path/to/denovo-map   
	cd /path/to/file
	ls -ltr     ##do files look like they copied over fine?? 
	rm /path/to/process-rads-directory/*.fq.gz    ##remove the files that ended up being copied over to the denovo-map)
	

#####3c. Looking at our demultiplexed dataset


Now that we've separated our dataset into our individuals, we can see how much data we have per individual. The fastest but not ideal way to look at this is simply by looking at the sizes of the .fq files that we now have for each individual. To get a list with file size information, type:

	ls -l

The files are huge!! and it's hard to read. So how can we look at them in a different unit (Megabyte? Terabyte?). We can find out by looking at the manual:

	man ls

(tip, scroll down to -h with spacebar; exit with 'q'). The you'll se that we can use that flag:

	ls -lh

and now we see more easily the sizes of our files. But can we order them? Let's pipe our list command with sort, using -rn (reverse numeric) to order from larger to smaller file:

	ls -l | sort -k 5 -rn

Uh oh.... looks like we have FIVE individuals with VERY little data!!! [Oh no!](http://giphy.com/gifs/movie-crying-johnny-depp-hAt4kMHnaVeNO). Well, it happens. So, we know we're likely to not use those individuals for our denovo-map step. But before we completely move on, if we want a more precise piece of information, we need read counts instead! First, gunzip your files:

	gunzip *.fq.gz

Then, count your reads (which will take a cuouple of minutes and will print to your screen):

	echo -e 'SAMPLE_ID_FULL\tNUM_READS'
	for file in ~/path/to/denovo-map/*.fq ##edit your path!!
	do
	echo -n $(basename $file .fq)$'\t'
	cat $file | grep '^@.*' | wc -l
	done

Here the syntax is more complicated, but in essence we are "grepping" (finding and grabbing) the first characters of each line of sequence for each file, and printing the numbers of lines with those characters that occur in each file... i.e., the number of reads for each individual.

We can see that we have those five individuals that have ZERO reads, and a couple of more that have very few reads.... However, for now let's only eliminate the individuals that have zero reads, to do denovo-map without those. 

	mkdir zero-reads
	mv Ch_327* zero-reads	
(and repeat that for all other four individuals). Check that you did this correctly by seeing there are no files with 0kb on your current directory and by ls -ls on the zero-reads directory to make sure all of them have zero reads. 


FYI: if you see from the file names, one of our individuals is named VERY differently from the others... that's because that individual sequence is a stray from another project. We don't want to run denovo-map with this individual, because it will skew our results. so, DELETE those files from the denovo-map directory:

	rm TNHC*

***NOTE:*** deleting files on the shell is much more dangerous than how you normally do this on finder on your computer... it permanently deletes them!! no going back... so always be sure and be careful. 

 OK, now we're ready to start mapping!!! 
 
	
###4. Running denovo map

#####4.1. Setting up your input file
You need to make a list (for example with ls -l) of all of your files because stacks needs as input the name of every single sequence in your dataset for running denovo_map.pl. In order to do this, we will copy/paste into a simple text editor all of your file names, and then run a chain of search/replace arguments to have this ready to input. 

First, keep only a single column with your file names. Then, add -s ./ to the begining of each line:

![find-replace-1](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/2-denovo-map/1-find-replace.png)

Using a series of finds that have to do with the first characters on each line:

![find-replace-2](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/2-denovo-map/2-find-replace.png)

Then, add a "line continues" to the end of each line (so that the shell knows to continue reading all lines of sequence):

![find-replace-3](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/2-denovo-map/3-find-replace.png)

And finally delete the last \ from the last line of sequence so that the program knows to STOP reading! (otherwise it will wait until you give it the last command...). The final file should look like [this](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/2-denovo-map/sequences-list.txt). 

***NOTE:*** specifically for the uniandes cluster, you need to specific the COMPLETE PATH every time you do anything (most clusters will find the sequence within the directory you're running the script from, simplifying paths). So, do a final search/replace to add entire path to file before each sequence. 



#####4.2. Setting up the denovo-map specifications to run a first time.

We can now add the specifications to denovo_map.pl BEFORE the list of sequences, and it should be like this:

	module load stacks/1.35

	mkdir ./path/to/denovo-map/denovo-test-1/

	denovo_map.pl -T 8 -m 2 -M 3 -n 2 -S -b 2 -o ./path/to/denovo-map/denovo-test-1/ \
	-s ./path/to/denovo-map/Ab_372.1.fq \
	-s ./path/to/denovo-map/Ab_372.2.fq \
Etc.... (as in, followed by all your other sequences). So, the final file should look like [this](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/2-denovo-map/denovo-map-test-2.txt), though you need to modify your personal information and path. 


#####4.3. Setting up denovo-map permutations for parameter testing

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


----------------------------------------------
----------------------------------------------

###5. Running *populations*

#####5.1. Setting up your popmap file
The popmap file is simply a list of all your individuals (i.e. the file names) similar to what you had to input for denovo_map, but without the file termination (as in, excluding fastq.gz). This list is in the first column, and the second column has the population that we've assigned to each of these individuals. In this case, we only have three populations, so we've named them 1, 2, and 3; though we could also name them as a short string of characters (such as ER, CH, AB). I've built the [popmap file](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/3-SNP-matrix/popmap-Stefania-c.txt) for you! But you should try to make this on your own for practice. 

#####5.2. Setting up your populations parameters and permutations
In the [populations](http://catchenlab.life.illinois.edu/stacks/comp/populations.php) program you essentially filter your SNP matrix and decide which outputs you want. 

You can start filtering your dataset here, though you can do much more filtering later on in something like R, sometimes you want to avoid dealing with a giant dataset of 200 thousand SNPs when you're sure you will filter out most of it due to missing data!! So we can set some baseline filters and use different parameter settings to see what our SNP output looks like. 

One thing we're sure we want for this dataset specifically is SNPs that are present in all three of our populations (if we had many more populations, and perhaps only few individuals or high divergence of some of them, then we might want to set this up differently). For this, we set the -p flag at 3. the other main parameter to vary is -r, which is the percent missing data you will allow per population to keep each SNP. In this case, we will start by setting a low bar of 40% (-r 0.4). The other important flag is --write_random_snp... since we don't have a reference genome, and we can have several snps within a read, we (for now) only want to keep a single snp per read so that we can "assume" our snps are not highly linked. However, we want those snps kept to be random (not the first or the last) in case something odd went on during Illumina sequencing (which has been observed....). The last few flags are the types of output we want. 


I have set up a [populations file](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/3-SNP-matrix/pops_r6_stef.sh) that's ready to go for p=3, r=0.6. You should try to set it up other settings. In this case, we will try out what happens when we set -r 0.4, -r 0.5, and -r 0.6. How does the SNP matrix look? too few SNPs with the highest threshold? too many with the lower one?? 





###6. Other important things to do for "real" datasets
There are MANY things you should do to make your your data is cleaned appropriately, before running actual population analyses, or at least to do as good a job as you can. Some of things include:

#####a. Running rxstacks for filtering possible paralogs
This [program](http://catchenlab.life.illinois.edu/stacks/comp/rxstacks.php) (in stacks)uses population haplotype information on a genome-wide level to infer which loci are having too many haplotypes given your dataset, and so filter out potential paralogs that were stacked up as orthologs. This re-runs stacks almost entirely, so it takes a while to run. 

#####b. Filtering Linked Loci and extremely rare alleles
Even though you don't have a reference genome, you should still filter out loci that seem to be linked. You do this after you have your SNP matrix and the linkage estimation is by SNP (allele) instead of using the entire sequence, which makes computation faster. You should also filter out minor allele frequencies of 0.05, which seems to be the standard for most of these SNP approaches. However, this runs the risk of overfiltering alleles that may be informative, for example, for inferring population expansion events. These filters are usually set in R, with different packages, for example in [genetics](http://svitsrv25.epfl.ch/R-doc/library/genetics/html/LD.html) or with [pegas](https://cran.r-project.org/web/packages/pegas/pegas.pdf).

#####c. Filter out weird Blast searches
It's quite likely your data is contaminated with things that are obviously not your species of interest. So you can do a massive blast search and filter out anything that is an obvious hit to something that is not even close (for example, bacteria if you're sequencing a frog!). 

###7. Running basic population analyses in R 
We will do some basic analyses in the R package [adegenet](http://adegenet.r-forge.r-project.org/), which is a fantastic package to get started, mostly because the developers have put a LOT of effor into making really handy tutorials and documentation.  

This, for example, is one of the most useful [basic tutorials](http://adegenet.r-forge.r-project.org/) I've ever found! 

When you try [this R script](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/4-quick-R-analyses/PCA_stef.R) for the first time with [this SNP matrix](https://github.com/pesalerno/Genomics-Intro-workshop/blob/master/4-quick-R-analyses/Stef_r7_1022-c.stru), it looks like some individuals may have too much missing data... can you identify them, remove them, remove those individuals, and rerun the script? do things look better/more logical? 

