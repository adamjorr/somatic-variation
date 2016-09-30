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
 5. Samtools version 1.2 and 1.3.1 or later (version 1.3 and above crashes Stampy, so a version of 1.2 must also be installed.)
 6. Stampy version 1.0.29 (r3755) or later.
 7. Bcftools version 1.3.1 or later
 8. HTSlib version 1.3.1 or later

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

Input File Requirements
-----------------------
For creation of a psuedo-reference:
 * A reference to a closely-related species in FASTA format
 * Paired sequencing reads in FASTQ format residing in a single directory, with reads from the first mate in a file labelled with an R1 and reads from the second mate in a file labelled with an R2. For example, reads_R1.fq and reads_R2.fq

Experimental Overview
---------------------
WIP (Explain the conditions surrounding the experiment here)


See the section on replicating our experiment for more information.


Pipeline Documentation
======================

Creating A Reference From a Similar Species with Iterative Mapping
------------------------------------------------------------------

Do kmer correction with:
```bash
bash clean_reads.sh ../data/
```

Then align to a reference with bwa:
```bash
bash bwa_aligner.sh ref.fa ./cleaned_reads/sliced/ bwa1.bam
```

Next, realign with Stampy:
```bash
bash stampy_realigner.sh ref.fa bwa1.bam stampy1.bam
```

Now obtain a consensus with:
```bash
bash create_consensus.sh -o consensus.fa ref.fa stampy1.bam
```

This process can be repeated, using the new consensus as a reference.

###clean_reads.sh
Usage: clean_reads.sh [-t THREADS] [-d DEST_DIRECTORY] [-k KMER_SIZE] [-f FILE_PATTERN] [-s SEARCH_STRING] [-r REPLACE_STRING] [-m MAX_MEMORY] [-c COVERAGE] READ_DIRECTORY

This script uses Rcorrector and Khmer to clean reads found in READ_DIRECTORY. The -m flag can be used to change the amount of memory given to Khmer. The -k flag changes the kmer size to use. The -f flag is a pattern that finds reads in READ_DIRECTORY matching that pattern; by default, it is `*R1*`. The -s and -r flags are the search and replace strings used to find the file containing the mates specified in the file found by the FILE_PATTERN. The default values are R1 and R2. For example, if read mates are specified only by 1 or 2 and not by R1 and R2, -s should be 1 and -r should be 2, so that the mates of reads1.fq are properly found in reads2.fq. The -c flag is an approximate coverage cutoff to filter on using Khmer and defaults to 40000.

###bwa_aligner.sh
Usage: bwa_aligner.sh [-t THREADS] [-p RG_PLATFORM] -r reference.fa -o out.bam data/

This script uses BWA to align reads in the given directory to the given reference and outputs a bamfile. The -p flag can be used to set the PLATFORM flag in the resultant BAM file; the default value is ILLUMINA.

###stampy_realigner.sh
Usage: stampy_realigner.sh [-t THREADS] [-d TMPDIR] [-q QUAL] reference.fasta data.bam out.bam

This script uses Stampy to realign unmapped reads and reads lower than QUAL in data.bam and outputs to out.bam. The given reference should be the same used to generate data.bam .
A temporary directory to be used can be specified with -d; the default is /tmp/.

Warning: Stampy tends to crash when Samtools 1.3 is used. Ensure that Samtools 1.2 is at the top of the path before running this script.

###create_consensus.sh
Usage: create_consensus.sh [-t THREADS] [-f FILTER] [-c CHAINFILE] [-b BCFTOOLSFILE] [-o OUTPUT] reference.fasta in.bam

This script uses bcftools and tabix to create a consensus in fasta format from a BAM file aligned to a reference. A filter is applied to the genotypes called by bcftools which can be changed with -f to give the desired behavior. A chain file that describes the differences between the reference and the output is generated and called consensus.chain by default. This can be changed with the -c flag. A file containing the genotype calls generated by BCFtools is generated and called bcftools_calls.vcf.gz by default; this can be changed with the -b flag.

Calling Somatic Variants and Filtering Results
----------------------------------------------
WIP

Experimental Replication
========================
WIP. Exact steps to replicate our analysis













