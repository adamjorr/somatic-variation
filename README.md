A Pipeline For the Detection of Somatic Variants
================================================
Herein we describe our pipeline for detecting somatic mutants. 

Many scripts in this pipeline use temporary files that are prone to clobbering, so don't attempt to run these scripts in parallel unless you do so in another working directory.

Software Requirements
---------------------
 1. [khmer](https://github.com/dib-lab/khmer) version 2.0+67.g136bb3d.dirty or later.
 2. [Rcorrector](https://github.com/mourisl/rcorrector) commit `eba3f79` or later.
 3. GNU parallel version 20160722 or later.
 4. BWA version 0.7.15-r1140 or later.
 5. Samtools version 1.3.1 or later.
 6. Stampy version 1.0.30 or later.
 7. Bcftools version 1.3.1 or later
 8. HTSlib version 1.3.1 or later
 9. GATK version 3.7 or later

Software Installation
---------------------
###khmer
```bash
pip install khmer
```

###Rcorrector
```bash
git clone https://github.com/mourisl/rcorrector.git
cd rcorrector
make
```

###GNU parallel
Many distributions have packages to install parallel. For Ubuntu:
```bash
sudo apt-get install parallel
```

###BWA
```bash
git clone https://github.com/lh3/bwa.git
cd bwa
make
```

###Samtools, Bcftools, and HTSlib
Visit the [htslib website](www.htslib.org/download/) for more information. In summary:

Download the release from www.htslib.org/download/, extract the archive, and enter the directory. Then,
```bash
make
make prefix=/where/to/install install
```

###Stampy
Stampy requires the user to register before downloading the software. Visit [this link](http://www.well.ox.ac.uk/software-download-registration) to do so. Once this is done, use the link provided in the registration email to download the archive. Extract it, enter the directory, then build the software with make.
```bash
make
```

###GATK
The GATK also requires registration before downloading the software. Visit [this link](https://software.broadinstitute.org/gatk/download/) and follow the instructions. Once you have accepted the license and downloaded the software archive, extract it and put the .jar file somewhere on your path.


Input File Requirements
-----------------------
For creation of a psuedo-reference:
 * A reference to a closely-related species in FASTA format
 * Paired sequencing reads in FASTQ format residing in a single directory, with reads from the first mate in a file labelled with an R1 and reads from the second mate in a file labelled with an R2. For example, reads_R1.fq and reads_R2.fq. Advanced users can get around this requirement using options to the `clean_reads.sh` script.

Experimental Overview
---------------------
TODO (Explain the conditions surrounding the experiment here)


See the section on replicating our experiment for more information.


Pipeline Summary
================

Creating A Reference From a Similar Species with Iterative Mapping
------------------------------------------------------------------

Do kmer correction with:
```bash
bash clean_reads.sh -i ../data/
```

Then align to a reference with bwa:
```bash
bash bwa_aligner.sh -r ref.fa -o bwa1.bam -i ./cleaned_reads/sliced/
```

Next, realign with Stampy:
```bash
bash stampy_realigner.sh -r ref.fa -i bwa1.bam -o stampy1.bam
```

Now obtain a consensus with:
```bash
bash create_consensus.sh -o consensus.fa -r ref.fa -i stampy1.bam
```

This process can be repeated, using the new consensus as a reference to converge on a reference specific to the individual sequenced.

Calling Somatic Variants 
------------------------
Variant calling can be performed using your favorite software. We provide a script that works with short Illumina reads and uses the GATK HaplotypeCaller to call variants.

To call variants:
```bash
bash gatkcaller.sh -r reference.fa -i data.bam -o variants.vcf
```

Filtering Results
------------------
To filter our VCF using any replication information, use `filt_with_replicates.pl`.

```bash
perl filt_with_replicates.pl <variants.vcf >filtered.vcf
```

Building A Tree
---------------

To turn this vcf into a FASTA file of all the variable sites, use `vcf2fa.sh`.
```bash
bash vcf2fa.sh <filtered.vcf >filtered.fa
```

To turn this FASTA file where each diploid genotype is one site into a FASTA file where each genotype is represented as 2 sites, use `diploidify.py`.
```bash
python2 diploidify.py <filtered.fa >doubled.fa
```

Now you can use your favorite program to create a tree using this FASTA file.

The filtering and tree construction can be simplified to one step if you have RAxML installed by using the `vcf2tree.sh` script.
```bash
bash vcf2tree.sh -i variants.vcf -o tree.nwk -g 3
```

We provide a simple R script to plot the tree and save it as a PDF. There are many tools that can visualize newick trees. This one will plot and label the tips using the first string of numbers that appears in its name.

```bash
Rscript plot_tree.R tree.nwk out.pdf
```


Experimental Replication
========================
```bash
do_experiment.sh folder_containing_reads/
```

Script Documentation
====================
##clean_reads.sh
Usage: clean_reads.sh [-t THREADS] [-d DEST_DIRECTORY] [-k KMER_SIZE] [-f FILE_PATTERN] [-1 FIRSTMATE] [-2 SECONDMATE] [-m MAX_MEMORY] [-c COVERAGE] -i READ_DIRECTORY

 * **-t:** number of threads to use [48]
 * **-d:** directory to put output in [./cleaned_reads]
 * **-k:** size of kmer to use [32]
 * **-f:** pattern matched to find reads ['*.fastq']
 * **-1:** pattern to find first mate in pair [R1]
 * **-2:** pattern to find second mate in pair [R2]
 * **-m:** max memory to give to khmer [64e9]
 * **-c:** max coverage to allow [40000]
 * **-i:** input directory to search for reads to clean

This script uses Rcorrector and Khmer to clean reads found in READ_DIRECTORY.
The -d option changes which directory to put output in and defaults to ./cleaned_reads .
The -k option changes the kmer size to use and defaults to 32.
The -f flag is a pattern that finds reads in READ_DIRECTORY matching that pattern; by default, it is '*.fastq' .
The -1 and -2 flags are the search and replace strings used to find the file containing the mates specified in the file found by the FILE_PATTERN. The default values are R1 and R2. For example, if read mates are specified only by 1 or 2 and not by R1 and R2, -s should be 1 and -r should be 2, so that the mates of reads1.fq are properly found in reads2.fq.
The -m flag controls how much memory is allocated to khmer and defaults to 64e9.
The -c flag is an approximate coverage cutoff to filter on using Khmer and defaults to 40000.
The -i option specifies the directory to search for reads to correct.

##bwa_aligner.sh
Usage: bwa_aligner.sh [-t THREADS] [-d TMPDIR] [-p RG_PLATFORM] [-q FILEPATTERN] [-1 FIRSTMATE] [-2 SECONDMATE] [-o out.bam] -r reference.fa -i data/

 * **-t:** number of threads to use [48]
 * **-d:** temporary directory to use
 * **-p:** PLATFORM flag for BAM file [ILLUMINA]
 * **-q:** pattern to find reads ['*.fastq']
 * **-1:** pattern to find first mate in pair [R1]
 * **-2:** pattern to find second mate in pair [R2]
 * **-o:** BAM file to write alignment to [STDOUT]
 * **-r:** reference to align to with BWA
 * **data/:** directory to search for reads

This script uses BWA to align reads in the given directory to the given reference and outputs a bamfile.
The -p flag can be used to set the PLATFORM flag in the resultant BAM file; the default value is ILLUMINA.
The -1, -2, and -q flags change how the script identifies the files containing reads to align.
-q is a pattern that is searched for to identify reads with the default value '*.fastq' .
To find files ending with '.fq', for example, use -q '*.fq'.
Remember to use quotes so the shell doesn't expand the wildcard.
The -1 and -2 options control how the script finds which of the FASTQ files given are the forward and reverse mate.
By default, -1 is R1 and -2 is R2.
The script will identify the first mates by looking for which FASTQ files contain the value of -1 in their names.
It will then match these files with their mates by substituting the value of -2.
The pairs of reads should be in the same directory and have the same name except for the substitution of the value of -1 for the value of -2.
For example, in a data directory with reads_R1.fastq and reads_R2.fastq, by default the script will identify the two files by their suffixes, identify the file containing R1 as the first set of reads and substitute R1 with R2 to identify the file containing the second set of reads.

##stampy_realigner.sh
Usage: stampy_realigner.sh [-t THREADS] [-d TMPDIR] [-o out.bam] -r reference.fasta -i data.bam

 * **-t:** number of threads to use [48]
 * **-d:** temporary directory to use
 * **-o:** BAM file to write alignment to [STDOUT]
 * **-r:** reference to align to with BWA
 * **-i:** input BAM file to remap

This script uses Stampy to realign poorly mapped reads (from Stampy's --bamkeepgoodreads option) and outputs to `out.bam`. The given reference should be the same used to generate `data.bam`.
A temporary directory to be used can be specified with -d.

##create_consensus.sh
Usage: create_consensus.sh [-t THREADS] [-f FILTER] [-c CHAINFILE] [-b BCFTOOLSFILE] [-o OUTPUT] reference.fasta in.bam

This script uses bcftools and tabix to create a consensus in fasta format from a BAM file aligned to a reference. A filter is applied to the genotypes called by bcftools which can be changed with -f to give the desired behavior. A chain file that describes the differences between the reference and the output is generated and called consensus.chain by default. This can be changed with the -c flag. A file containing the genotype calls generated by BCFtools is generated and called bcftools_calls.vcf.gz by default; this can be changed with the -b flag.













