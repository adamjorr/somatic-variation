A Pipeline For the Detection of Somatic Variants
================================================
Herein we describe our pipeline for detecting somatic mutants. 

Many scripts in this pipeline use temporary files that are prone to clobbering, so don't attempt to run these scripts in parallel unless you do so in another working directory.

Overview
-------
```bash
clean_reads.sh -i /dir/with/reads #correct reads with Rcorrector, then use khmer to filter out reads with repetitive kmers.
ngm_aligner.sh -r original_ref.fa -i ./cleaned_reads/sliced/ -o original_align.bam #align the reads in the input directory with ngm
create_consensus.sh -r original_ref.fa -i original_align.bam -o better_ref.fa #create a consensus sequence to use as an individualized reference
ngm_aligner.sh -r better_ref.fa -i /dir/with/reads -o better_alignment.bam #align uncorrected reads to the individualized reference
gatkcaller.sh -r better_ref.fa -i better_alignment.bam -o variants.vcf #use gatk to call variants
bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <=40' variants.vcf | bcftools view -m2 -M2 -v snps -o filtered.vcf #filter variants and grab only snps
vcf2tree.sh -i filtered.vcf -o tree.nwk #concatenate the snps the generate a tree
```

Important Filters
-----------------
These scripts have important filters. Their behavior is summarized here.
 * `create_consensus.sh` - By default, ignore sites where number of alt alleles equals 1 while creating consensus sequence.
   - **-f** can be used to change this with any valid bcftools filter expression
 * `filt_with_replicates.pl` - Samples with discordant replicates have their genotypes changed to ./.
   - **-s:** no missing genotypes allowed (only emit sites where every sample's replicates agree)
   - **-m:** change non-matching genotypes to the one that matches the majority
 * `diploidify.py`
   - removes sites that require two mutations (like C/C -> T/T)
   - removes sites that have 3 or more alleles
   - **-v** outputs only variable sites

We also recommend filtering your final variants with a depth and heterozygosity filter, as well as ignoring variants near indels (the `-g` option).
These can be symptoms of alignment errors.
We also remove any non-snps. This is done with the following command:
```bash
bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <=40' repfiltered.vcf | bcftools view -m2 -M2 -v snps -o filtered.vcf
```

Software Requirements
---------------------
 1. [khmer](https://github.com/dib-lab/khmer) version 2.0+67.g136bb3d.dirty or later to use `clean_reads.sh`.
 2. [Rcorrector](https://github.com/mourisl/rcorrector) commit `eba3f79` or later to use `clean_reads.sh`.
 3. [NextGenMap](http://cibiv.github.io/NextGenMap/) version 0.5.2 or later.
 4. GNU parallel version 20160722 or later.
 5. Samtools, Bcftools, and HTSlib version 1.3.1 or later.
 6. GATK version 3.7 or later to use `gatkcaller.sh`
 7. RAxML version 8.0.0 or later to use `vcf2tree.sh`
To fully replicate our experiment, you will also need:
 8. Bedtools
 9. liftOver
 10. VCFtools

Notes For Installing Stuff with Conda
--------------------------------------
```bash
conda install -n somatic-variation -c bioconda/label/gcc7 khmer rcorrector nextgenmap samtools bcftools gatk4 raxml bedtools ucsc-liftover vcftools perl-vcftools-vcf parallel whatshap perl-bio-cigar picard pyvcf biopython r-ape r-phytools r-tidyverse r-rcolorbrewer r-vcfr r-phangorn r-doparallel r-dfoptim intermine bioconductor-karyoploter bioconductor-variantannotation bioconductor-iranges bioconductor-genomicranges bioconductor-genomicfeatures bioconductor-rtracklayer pysam pybedtools 
```

Software Installation
---------------------
### khmer
```bash
pip install khmer
```

If you're having trouble with pip or prefer conda, you can also install khmer with:

```bash
conda install khmer
```


### Rcorrector
```bash
git clone https://github.com/mourisl/rcorrector.git
cd rcorrector
make
```

### NextGenMap
Ensure your conda channels are properly configured
```bash
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
```
Then install
```bash
conda install nextgenmap
```

### GNU parallel
Many distributions have packages to install parallel. For Ubuntu:
```bash
sudo apt-get install parallel
```

### Samtools, Bcftools, and HTSlib
Visit the [htslib website](www.htslib.org/download/) for more information. In summary:

Download the release from www.htslib.org/download/, extract the archive, and enter the directory. Then,
```bash
make
make prefix=/where/to/install install
```

### GATK
The GATK requires registration before downloading the software. Visit [this link](https://software.broadinstitute.org/gatk/download/) and follow the instructions. Once you have accepted the license and downloaded the software archive, extract it and put the .jar file somewhere on your path.

### RAxML
Detailed installation instructions are available on [github](https://github.com/stamatak/standard-RAxML). To summarize, run:
```bash
git clone https://github.com/stamatak/standard-RAxML
cd standard-RAxML
make -f Makefile.PTHREADS.gcc
rm *.o
```

A different compilation configuration may allow RAxML to run faster based on your computer's configuration, so be sure to look at the repository for more information.

### BEDTools
Download the latest release from [github](https://github.com/arq5x/bedtools2/), extract it, enter the directory and `make`.
```bash
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2
make
```
You should then add the binaries to your path.

### liftOver
Precompiled binaries are available from UCSC. To download liftOver,
```bash
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
```
Then add the binary to your bath.

### VCFtools
Visit the VCFtools [website](https://vcftools.github.io/examples.html) and download the archive of the latest release. Extract the directory, set your PERL5LIB evironment variable to contain the path to Vcf.pm, then compile the program with `./configure`, `make`, and `make install`. To summarize:
```bash
tar -zxvf vcfools.0.X.XX.tar.gz
export PERL5LIB=/path/to/your/vcftools-directory/src/perl/
cd vcftools/ 
./configure	
make 
make install
```

Input File Requirements
-----------------------
For creation of a psuedo-reference:
 * A reference to a closely-related species in FASTA format
 * Paired sequencing reads in FASTQ format residing in a single directory, with reads from the first mate in a file labelled with an R1 and reads from the second mate in a file labelled with an R2. For example, reads_R1.fq and reads_R2.fq. This requirement can be circumvented using options to the `clean_reads.sh` script.

To filter variants using replicated data, you must have replicate samples. These samples should be named such that the replicate identifier is the last character in the string. For example, sample M1a, M1b, and M1c are all replicates of sample M1. M2a, M2b, and M2c are all replicates of M2, which is a different sample than M1. This requirement can be circumvented by ordering the samples in the VCF such that replicates are next to one another and using the `-g` argument to `filt_with_replicates.pl`. This can only be done if you have an equal number of replicates for all your samples.

Our Use Case
------------
We investigate the following specific questions:
 * Can we reconstruct the topology of a tree using only genetic data sampled from its leaves?
 * What is the somatic mutation rate of a long-lived eucalyptus?

To answer these questions, we sampled in triplicate 8 branches of a eucalyptus tree.
This tree does not have a reference genome, which is why we attempt to clean up the genome before aligning and calling variants.
If you have a high-quality reference genome for your data, feel free to skip that step.

See the section on replicating our experiment for more information.

Pipeline Summary
================

Creating A Reference From a Similar Species with Iterative Mapping
------------------------------------------------------------------

We first clean the reads to avoid fitting errors.

Do read correction with:
```bash
bash clean_reads.sh -i /dir/with/reads
```

Then align to a reference with NextGenMap:
```bash
ngm_aligner.sh -r original_ref.fa -i ./cleaned_reads/sliced/ -o original_align.bam
```

Now obtain a consensus with:
```bash
create_consensus.sh -r original_ref.fa -i original_align.bam -o better_ref.fa
```

This process can be repeated, using the new consensus as a reference to converge on a reference specific to the individual sequenced.

To avoid interfering with any assumptions made by the variant caller, we recommend aligning the uncorrected reads before calling variants.
```bash
ngm_aligner.sh -r better_ref.fa -i /dir/with/reads -o better_alignment.bam
```



Calling Somatic Variants 
------------------------
Variant calling can be performed using your favorite software. We provide a script that works with short Illumina reads and uses the GATK HaplotypeCaller to call variants.

To call variants:
```bash
bash gatkcaller.sh -r better_ref.fa -i better_alignment.bam -o variants.vcf
```

Filtering Results
------------------
To filter our VCF using available replication information, use `filt_with_replicates.pl`.
By default this script assumes replicates of a sample will be identically named except for a replicate identifier as the last character.
You can use the `-g` option to tell the script that you have your VCF file ordered such that the replicates are next to each other.
For example, `-g 3` will assume that the first 3 samples are replicates of one another, and every group of 3 samples after that are replicates of each other, and so on.
```bash
perl filt_with_replicates.pl <variants.vcf >repfiltered.vcf
```

We recommend following this with a depth and heterozygosity filter and removing any non-snps. This can be accomplished by running
```bash
bcftools filter -g 50 -i 'DP <= 500 && ExcessHet <=40' repfiltered.vcf | bcftools view -m2 -M2 -v snps -o filtered.vcf
```
This gives us a filtered VCF of our 'best' variant calls. 

Building A Tree
---------------

To turn this vcf into a FASTA file that concatenates all the variable sites, use `vcf2fa.sh`.
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
bash vcf2tree.sh -i variants.vcf -o tree.nwk
```

We provide a simple R script to plot the tree and save it as a PDF. There are many tools that can visualize newick trees. This one will plot and label the tips using the first string of numbers that appears in its name.

```bash
Rscript plot_tree.R tree.nwk out.pdf
```

Experimental Replication
========================

Extra Filters
-------------
We find the quality of our variants increases when filtering out variants in repetitive regions.
To do so, we liftOver repetitive region annotations to our individualized reference and filter with bcftools.
For example,

```bash
liftOver -gff mask.gff better_ref.fa.chain better_ref_mask.gff /dev/null
```

We can then convert this to a BED file and exclude any positions in the VCF that overlap a site in the BED file with VCFtools.
We use the chainfile automatically generated by `create_consensus.sh`.
```bash
gff2bed <better_ref_mask.gff >better_ref_mask.bed
vcftools --vcf filtered.vcf --exclude-bed better_ref_mask.bed --remove-filtered-all --recode --recode-INFO-all --stdout >filtered_bed_excluded.vcf
```

Makefile
--------
We provide a Makefile to allow for replication of our experiment and as an example of a pipeline.
We provide default variables, but they can be overwritten by setting the variable while calling `make`; for example, `make REF=myreference.fa`.

```bash
make
```

Script Documentation
====================
## clean_reads.sh
Usage: clean_reads.sh [-t THREADS] [-s SLICE_THREADS] [-d DEST_DIRECTORY] [-k KMER_SIZE] [-f FILE_SUFFIX] [-1 FIRSTMATE] [-2 SECONDMATE] [-m MAX_MEMORY] [-c COVERAGE] -i READ_DIRECTORY

 * **-t:** number of threads to use [48]
 * **-s:** threads to use during slicing [THREADS/6]
 * **-d:** directory to put output in [./cleaned_reads]
 * **-k:** size of kmer to use [32]
 * **-f:** pattern matched to find reads [.fastq]
 * **-1:** pattern to find first mate in pair [R1]
 * **-2:** pattern to find second mate in pair [R2]
 * **-m:** max memory to give to khmer [64e9]
 * **-c:** max coverage to allow [40000]
 * **-i:** input directory to search for reads to clean

This script uses Rcorrector and Khmer to clean reads found in READ_DIRECTORY. Rcorrector is used to remove errors and Khmer is used to remove excessively repetitive reads.
The -d option changes which directory to put output in and defaults to ./cleaned_reads .
The -k option changes the kmer size to use and defaults to 32.
The -f flag is a pattern that finds reads in READ_DIRECTORY matching that pattern; by default, it is '\*.fastq' .
The -1 and -2 flags are the search and replace strings used to find the file containing the mates specified in the file found by the FILE_PATTERN. The default values are R1 and R2. For example, if read mates are specified only by 1 or 2 and not by R1 and R2, -s should be 1 and -r should be 2, so that the mates of reads1.fq are properly found in reads2.fq.
The -m flag controls how much memory is allocated to khmer and defaults to 64e9.
The -c flag is an approximate coverage cutoff to filter on using Khmer and defaults to 40000.
The -i option specifies the directory to search for reads to correct.
The -s option specifies the number of threads to use during khmer slicing. This will use MAX_MEMORY x SLICE_THREADS amount of memory, so make sure your system supports at least that much memory or the program will crash.

### What it does
Uses Rcorrector and Khmer to clean reads found in READ_DIRECTORY. Rcorrector is used to remove errors and Khmer is used to remove excessively repetitive reads.

### Why do we use it
Later in the pipeline, we want to estimate a genome using the reads we've generated and the reference of a closely related species.
To reduce the chance of fitting errors during the creation of a new reference, we remove sequencing errors first with Rcorrector.
We then remove extremely high-coverage reads with Khmer, as these are likely erroneous and may cause errors if they map to the wrong place in the genome.

### Expected output
A file containing the khmer graph, called khmer_count.graph.
Two folders: 
 * **corrected/** containing the corrected read pairs.
 * **sliced/** containing the read pairs after filtration with Khmer, as well as reads which passed filtering but their mate did not.

Output reads will have the same name as the input reads, but corrected reads will have a .cor.fq suffix and sliced reads will have a .cor_sliced.${FILE_SUFFIX} suffix.

## slice-paired-reads-by-coverage.py
Usage: slice-paired-reads-by-coverage.py [-m MIN_COVERAGE] [-M MAX_COVERAGE] INPUT_GRAPH INPUT_READS1 INPUT_READS2 OUTPUT_READS1 OUTPUT_READS2 [OUTPUT_SINGLETONS]

 * **-m:** remove reads with approximate coverage below this value [none]
 * **-M:** remove reads with approximate coverage above this value [none]
 * **INPUT_GRAPH:** khmer graph to use
 * **INPUT_READS1:** first end of paired input reads
 * **INPUT_READS2:** second end of paired input reads
 * **OUTPUT_READS1:** first end of paired output reads
 * **OUTPUT_READS2:** second end of paired output reads
 * **OUTPUT_SINGLETONS:** file to ouput single reads to if their mate fails filtering [none]

This is a modified version of khmer's [slice-reads-by-coverage.py](https://github.com/dib-lab/khmer/blob/master/sandbox/slice-reads-by-coverage.py) script. It can handle paired-end reads. Check out [this blog post](http://ivory.idyll.org/blog/2014-slice-reads-by-coverage.html) for more details on how the script works. In summary, it will remove reads which have a high or low expected coverage level. This can be useful for removing reads which are likely to be errors or reads likely to map to repetitive regions. All positional arguments are required except the singleton output file. One of -m or -M must be specified.

### What it does
Removes reads which have a high or low expected coverage level using khmer.

### Why do we use it
This is a modified version of a script included in the khmer package. The included version of the script does not keep any read pairing information. Since we have
paired reads and would like to preserve mate pair information, but also want to filter out highly repetitive reads (as described above), we use this script.

### Expected output
3 files: OUTPUT_READS1, OUTPUT_READS2, and OUTPUT_SINGLETONS, which are the first, second, and singleton reads. The singleton reads are those that passed filtering, but their mates did not.

## ngm_aligner.sh
Usage: bwa_aligner.sh [-t THREADS] [-d TMPDIR] [-p RG_PLATFORM] [-q FILEPATTERN] [-1 FIRSTMATE] [-2 SECONDMATE] [-o out.bam] [-s SENSITIVITY] -r reference.fa -i data/

 * **-t:** number of threads to use [48]
 * **-d:** temporary directory to use
 * **-p:** PLATFORM flag for BAM file [ILLUMINA]
 * **-q:** pattern to find reads ['\*.fastq']
 * **-1:** pattern to find first mate in pair [R1]
 * **-2:** pattern to find second mate in pair [R2]
 * **-o:** BAM file to write alignment to [STDOUT]
 * **-s:** SENSITIVITY parameter to pass to ngm [estimated from data]
 * **-r:** reference to align to with ngm
 * **data/:** directory to search for reads

This script uses ngm to align reads in the given directory to the given reference and outputs a bamfile.
The -p flag can be used to set the PLATFORM flag in the resultant BAM file; the default value is ILLUMINA.
The -1, -2, and -q flags change how the script identifies the files containing reads to align.
-q is a pattern that is searched for to identify reads with the default value '\*.fastq' .
To find files ending with '.fq', for example, use -q '\*.fq'.
Remember to use quotes so the shell doesn't expand the wildcard.
The -1 and -2 options control how the script finds which of the FASTQ files given are the forward and reverse mate.
By default, -1 is R1 and -2 is R2.
The script will identify the first mates by looking for which FASTQ files contain the value of -1 in their names.
It will then match these files with their mates by substituting the value of -2.
The pairs of reads should be in the same directory and have the same name except for the substitution of the value of -1 for the value of -2.
For example, in a data directory with reads_R1.fastq and reads_R2.fastq, by default the script will identify the two files by their suffixes, identify the file containing R1 as the first set of reads and substitute R1 with R2 to identify the file containing the second set of reads.

### What it does
Uses ngm to align reads in the given directory to the given reference and outputs a bamfile.

### Why do we use it
ngm is a fast and accurate way to map reads to a reference. It also calculates a sensitivity parameter that can be helpful when mapping to a divergent reference.

### Expected output
A bam file, specified by -o.

## create_consensus.sh
Usage: create_consensus.sh [-t THREADS] [-d TMPDIR] [-f FILTER] [-c CHAINFILE] [-b BCFTOOLSFILE] -o out.fa -r reference.fasta -i in.bam

 * **-t:** number of threads to use [48]
 * **-d:** temporary directory to use
 * **-f:** filter to apply for creating the consensus [AC==1]
 * **-c:** name of chainfile (for use with the liftover program) to generate [$outfile.chain]
 * **-b:** name of gzipped VCF to use to save bcftools' variant calls to [none]
 * **-o:** output file to write consensus FASTA.
 * **-r:** reference used to create input BAM file.
 * **-i:** input BAM file to create consensus from.

This script uses bcftools and tabix to create a consensus in fasta format from a BAM file aligned to a reference. A filter is applied to the genotypes called by bcftools which can be changed with -f to give the desired behavior. The default is to ignore sites where the number of alt alleles equals 1. A chain file that describes the differences between the reference and the output is generated with the same name as the output fasta file with an added .chain suffix by default. This can be changed with the -c flag. A file containing the genotype calls made by BCFtools is generated and can be saved for use with the -b option.

### What it does
Uses bcftools and tabix to create a consensus sequence in fasta format from a BAM file.

### Why do we use it
Since we don't have a proper reference for our data, we call a consensus from our alignment to generate a new reference that is hopefully closer to
the true reference for our species. Using this reference will improve the quality of an alignment produced using reads from our species.
A higher quality alignment will lead to higher quality variant calls.

### Expected output
 * A consensus reference specified by -o.
 * A chainfile describing the rearrangements between the input reference and the output consensus file. By default, this has the same name as the output consensus but with an additional .chain suffix.
 * When -b is used, a vcf file of variant calls used to generate the consensus.

## gatkcaller.sh
Usage: gatkcaller.sh [-t THREADS] [-p PICARD_CMD] [-d TMPDIR] [-g GATK_PATH] [-b BEDFILE] [-o OUTFILE] -r ref.fa -i data.bam

 * **-t:** number of threads to use [48]
 * **-p:** command to invoke picard on this system [picard]
 * **-d:** temporary directory to use
 * **-g:** path to GATK jar file [~/bin/GenomeAnalysisTK.jar]
 * **-b:** BED file of regions to exclude from variant calling [none]
 * **-o:** output file to write variants to [stdout]
 * **-r:** reference used to create input BAM file.
 * **-i:** input BAM file to call variants from.

A script implementing the [GATK best practices](https://software.broadinstitute.org/gatk/best-practices/) to call variants from the input BAM file. A BED file can be specified to exclude the regions described in it from variant calling, if desired.

### What it does
Calls variants from a BAM file using the GATK best practices.

### Why do we use it
Discovering variants from an alignment of reads is not trivial, so we use the GATK best practices pipeline to do so. Though we are interested in somatic variants,
we don't use the GATK's somatic variant caller because it is intended for detection of variation in human cancer.

### Expected output
A vcf file specified by -o.

## filt_with_replicates.pl
Usage: filt_with_replicates.pl [-s | --strict] [-m | --majority] [-g | --grouped int] [-h | -? | --help] <vcfile.vcf >filtered.vcf

 * **-s:** do strict filtering
 * **-m:** do majority-rule filtering
 * **-g:** replicates are found in groups of INT

This script filters a vcf using replicate data from the VCF. By default, it will consider all samples to be replicates of one another that have the same sample name except for the last character. Using -g, one can specify that the replicates are grouped together in the VCF. For example, -g 3 will assume that the first 3 samples are replicates of one another, the 2nd 3 samples are replicates of one another, and so on. A group of replicates will not pass the filter unless the genotypes of all replicates match. The -m flag relaxes this stipulation, and allows a group to pass filtering if the majority of replicates match. In this case, all genotypes of the replicate will be changed to the majority. If a replicate fails to pass filtering, all samples of that replicate are changed to have a ./. genotype. The -s flag alters filtering behavior such that a site is only emitted if all replicates pass filtering at that site.

### What it does
Filter a vcf using replicate information.

### Why do we use it
The GATK pipeline doesn't use information about replicates to determine the quality of variants. We have replicates and can use that information to reduce the
number of false positives with this script.

### Expected output
A filtered vcf file printed to stdout.

## vcf2fa.sh
Usage: vcf2fa.sh [-d tmp_dir/] [-i file.vcf] [-o out.fa]

 * **-d:** temporary directory to use [/tmp/ or equivalent]
 * **-i:** input vcf [stdin]
 * **-o:** output fasta file [stdout]

This script takes sites in the input VCF and outputs them to a fasta alignment of all the sites using IUPAC ambiguity for heterozygous sites.
Note that any downstream analysis should be done using ascertainment bias correction if possible, unless you know exactly what you're doing.

### What it does
Takes sites from a VCF file and generates a fasta alignment of these sites.

### Why do we use it
The easiest way to generate a phylogeny with any arbitrary software package is to give it a fasta-formatted alignment of sites.
We use this script to generate that alignment from our file of variants.

### Expected output
A fasta alignment specified by -o. If -o is not used, prints to stdout.

## diploidify.py
Usage: diploidify.py [-i INFILE] [-t IN_TYPE] [-o OUTFILE] [-p OUTTYPE] [-v] [-s]

 * **-i:** input alignment [stdin]
 * **-t:** input alignment type [fasta]
 * **-o:** output alignment [stdout]
 * **-p:** output alignment type [fasta]
 * **-v:** use to output only variable sites
 * **-s:** use to keep single-letter notation but keep other checks

A script to turn an alignment where each diploid genotype is represented as one site with IUPAC ambiguity codes into an alignment where each genotype is represented as 2 sites. This script also removes sites that require two mutations (like C/C -> T/T). It also removes sites that have 3 alleles. Using -s will keep these checks (ie remove triallelic sites and sites that require multiple mutations) but won't expand each site into two other sites, keeping the original notation. Using -v will output only variable sites.

### What it does
Takes a fasta alignment and doubles the number of sites, eliminating IUPAC ambiguity codes. This is useful if the ambiguity codes represent a diploid genotype and not actual ambiguity.

### Why do we use it
We have a diploid organism, and at each site the genotype of our diploid is represented with an IUPAC ambiguity code. However, most tree inferring software
will assume (correctly) that these codes represent ambiguity at a site. Thus, we use this script to remove ambiguiuty codes while retaining information about
mutation and differences between samples at any site.

### Expected output
An alignment specified by -o. If -o is not used, prints to stdout.

## vcf2tree.sh
Usage: vcf2tree.sh [-t THREADS] [-r raxmlHPC] [-i in.vcf] [-o out.nwk] -g GROUPBY

 * **-t:** threads to use for tree construction [12]
 * **-r:** how to invoke RAxML on this system [raxmlHPC]
 * **-i:** input vcf [stdin]
 * **-o:** output newick tree [stdout]
 * **-g:** replicates are in groups of this. See documentation for `filt_with_replicates.pl` for more information.

Takes a vcf file and filters it using the strict mode of `filt_with_replicates.pl`, creates a fasta alignment with `vcf2fa.sh`, then creates a diploidified version of the alignment with `diploidify.py`, and finally creates a newick tree using RAxML. -g is the same argument used in `filt_with_replicates.pl`. Use -r to tell the script how to invoke RAxML on this system.

### What it does
Creates a newick tree from a VCF file.

### Why do we use it
To generate a tree with RAxML.

### Expected output
A newick tree specified by -o. If -o is not specified, prints the tree to stdout.

## plot_tree.R
Usage: Rscript plot_tree.R tree.nwk plot.pdf

 * **tree.nwk:** name of tree file to plot
 * **plot.pdf:** name of pdf file to write plot to

R script to plot the given newick tree. The script will label the tips using the first string of numbers that appears in its name. Requires the ape and phytools R packages.

### What it does
Plots a newick tree.

### Why do we use it
We want to see a graphical representation of what a tree looks like.

### Expected output
A pdf file, given as the second argument to the script.

