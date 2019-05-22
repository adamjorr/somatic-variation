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

Notes For Installing with Conda
-------------------------------
```bash
conda install -n somatic-variation -c bioconda/label/gcc7 khmer rcorrector nextgenmap samtools bcftools gatk4 raxml bedtools ucsc-liftover vcftools perl-vcftools-vcf parallel whatshap perl-bio-cigar picard pyvcf biopython r-ape r-phytools r-tidyverse r-rcolorbrewer r-vcfr r-phangorn r-doparallel r-dfoptim intermine bioconductor-karyoploter bioconductor-variantannotation bioconductor-iranges bioconductor-genomicranges bioconductor-genomicfeatures bioconductor-rtracklayer pysam pybedtools 
```

Add to path in the environment, remove from path upon deactivation:
```bash
mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d
mkdir -p ${CONDA_PREFIX}/etc/conda/deactivate.d
cat >${CONDA_PREFIX}/etc/conda/activate.d/somatic-variation-modpath.sh <<EOF
export SOMATICVAROLDPATH=${PATH}
export PATH="$(readlink -f ./scripts/)":${PATH}
EOF

cat >${CONDA_PREFIX}/etc/conda/deactivate.d/somatic-variation-modpath.sh <<EOF
export PATH=${SOMATICVAROLDPATH}
unset SOMATICVAROLDPATH
EOF
```

Potential Problems
------------------
 * lines 71-79 of analysis/dng/filter/Makefile references some scripts to make figures. They need to be added to a subdirectory.

Doing the Analysis with the provided Makefiles
----------------------------------------------
The recommended way to run the analysis is with the makefiles provided in the
:open_file_folder: analysis/ and :open_file_folder: data/ directories. We provide
some helper scripts in the :open_file_folder: script/ directory that are sometimes
called by the Makefiles.


Complete File Structure after everything is made
------------------------------------------------
 * :open_file_folder: data/
 	* :open_file_folder: e_grandis/ - resources from the *E. grandis* genome
 		* ref.fa - *E. grandis* reference
 		* egrandis_201_repeatmask.gff3 - *E. grandis* repeat mask
 		* egrandis_201_genes.gff3 - *E. grandis* gene annotations
 	* :open_file_folder: raw/ - the raw reads
 	* :open_file_folder: mutated/ - reads with simulated mutations to determine the false negative rate
 * :open_file_folder: analysis/
 	* :open_file_folder: cleaned_reads/ - a folder containing high confidence reads for cleaning up reference
 	* :open_file_folder: e_mel_1/ - first round of reference clean up
 		* e_mel_1.bam - cleaned reads aligned to grandis reference
 		* e_mel_1.fa - consensus called from alignment
 	* :open_file_folder: e_mel_2/ - second round of reference clean up
 		* e_mel_2.bam - cleaned reads aligned to e_mel_1.fa
 		* e_mel_2.fa - consensus called from alignment
 	* :open_file_folder: e_mel_3/ - third round of reference clean up
 		* e_mel_3.bam - cleaned reads aligned to e_mel_2.fa
 		* e_mel_3.fa - consensus called from alignment
 		* e_mel_3.bed - e_mel_3 chromosome sizes in BED format
 		* e_mel_3_norepeats.bed - e_mel_3 chromosome sizes minus the repeatmask where regions are length >= 200.
 		* e_mel_3_norepeats.regions - e_mel_3_norepeats.bed reformatted to the 1-based inclusive format (chr:start-end)
 	* :open_file_folder: liftover/ - grandis annotations lifted over to coordinates of new references
 		* e_mel_1_repeatmask.gff3
 		* e_mel_2_repeatmask.gff3
 		* e_mel_3_repeatmask.gff3
 		* e_mel_3_repeatmask.bed
 		* e_mel_1_genes.gff3
 		* e_mel_2_genes.gff3
 		* e_mel_3_genes.gff3
 		* e_mel_3_genes.bed
 	* alignment.bam - raw sequence data aligned to the e_mel_3.fa reference
 	* alignment.dedup.bam - alignment.bam deduped with picard
 	* alignment.firstcalls.vcf - conservative calls made with Haplotypecaller as a first pass
 	* alignment.recal.table - recalibration table for BQSR
 	* alignment.recal.bam - deduped and recalibrated alignment
 	* e_mel_3_cds.bed - CDSs in e_mel_3 coordinates. Used for finding degenerate sites
 	* :open_file_folder: degenerate/ - finding degenerate sites to estimate model parameters
 		* e_mel_3_degenerate.sites
 		* e_mel_3_degenerate.bed
 	* alignment_norepeats.ad - DeNovoGear pileup of sites at e_mel_3_norepeats.regions
 	* :open_file_folder: dng
 		* :open_file_folder: gatk/ - initial call set to improve
 			* *scaffold_#*.haplotypecaller.vcf.gz - initial calls for the specified scaffold
 		* :open_file_folder: call/ - DeNovoGear calls
 			* *scaffold_#*.call.vcf.gz - calls from DeNovoGear
 		* allgatk.vcf.gz - concatenated GATK initial calls
 		* allcalls.vcf.gz - concatenated DeNovoGear calls
 		* denovos.vcf.gz - only denovo mutations from the DeNovoGear calls
 		* :open_file_folder: filter/ - further filters on the DeNovoGear callset
 			* genome.txt - genome description for BEDTools
 			* deduped.bed - variant sites that do not have another variant within 1000 nucleotides
 			* deduped.vcf.gz - vcf of variant sites that don't have another variant within 1000 nucleotides
 			* deduped.slop.bed - 500 bases on each side of each variant in deduped.bed
 			* bam/deduped_regions.bam - BAM of reads that mapped in deduped.slop.bed
 			* hplen.tsv.gz - tsv containing length of haplotype blocks each variant is in
 			* goodsites.vcf - vcf containing denovo sites in haplotype blocks over 500bp long
 			* goodsites.slop.bed - 100 bases on each side of the variants in goodsites
 			* goodsites.context.txt - list of base before, ref base, and base after each variant in goodsites
 		* :open_file_folder: optim/ - files related to fitting model parameters
			* *scaffold_#*.bed - degenerate sites on the specified scaffold
			* *scaffold_#*.bcftools.vcf.gz - bcftools called on the degenerate sites
			* loglike.jobs - file containing commands to run to obtain the log likelihood of the model
			* optimized_params.txt - file describing maximum likelihood parameters
 		* :open_file_folder: phasing/ - files related to phasing variants
			* :open_file_folder: input/ - input files for phasing
				* *scaffold_#*.vcf.gz - input variants to phase
			* :open_file_folder: phased/ - output phased files
				* *scaffold_#*.vcf.gz - whatshap phased variants
			* phased.vcf.gz - combined phased variants
 	* :open_file_folder: false_negative_rate/ - files for finding the false negative rate. This Makefile makes the mutated reads.
	 	* sites_to_mutate.bed - random sites to mutate
		* sites_to_mutate.vcf.gz - random sites genotyped using parameters estimated from bcftools calls and a star tree
		* sites_to_mutate.txt - genotyped random sites in a tsv format suitable for the mutation simulation script
		* :open_file_folder: mut_files/ - lists of mutations to create in each sample
			* mutations.tab - locations of mutations
			* mutfile.txt - master list of mutations
			* mut_M#.txt - sample-specific mutations
		* :open_file_folder: sed_scripts/*M##*/ - sed scripts for mutating the sample
			* mut_M*#*.txt.R1.sed - sed script for mutating the R1 reads
			* mut_M*#*.txt.R2.sed - sed script for mutating the R2 reads
		* alignment.bam - alignment of mutated reads to e_mel_3 reference
		* alignment.dedup.bam - alignment deduped with picard
		* alignment.firstcalls.vcf - first calls used to do BQSR
		* alignment.recal.table - recalibration table used to do BQSR
		* alignment.recal.bam - alignment after deduplication and BQSR
		* :open_file_folder: dng/ - finding the false negative rate for the denovogear pipeline
			* denovos.vcf.gz - denovo mutations called from the mutated sites
			* inrepeats.txt - mutation locations that are in the repeatmask
			* results.txt - tab delimited file showing each simulated mutation and whether or not it was detected
			* :open_file_folder: filter/ - same as the files in analysis/dng/filter but created from the alignment containing the mutated reads
			* :open_file_folderL phasing/ - same as the files in analysis/dng/phasing but created from the alignment containing the mutated reads
 	* :open_file_folder: false_positive_rate/ - files for estimating the false positive rate of the pipeline
	 	* :open_file_folder: trees/ - randomly generated trees to use to call variants
		* :open_file_folder: filter/ - removing clustered variants for each callset
		* :open_file_folder: goodsites/ - further filtering based on haplotype block length
		* :open_file_folder: denovos/ - denovos called for each randomly generated tree
		* allgatk.vcf.gz - variants previously called with GATK
		* num_false_pos.txt - number of variants called in each final call set
		* allgoodsites.bed - BED file containing every called site in all simulations
		* num_false_pos_no_originalcalls.txt - num_false_pos.txt but excluding any site that was called as variable using the true tree
		* originalcalls_overlaps.vcf.gz - original calls that overlap a site detected in any of the simulated calls
 	* :open_file_folder: gatk4/ - calling variants with GATK to determine whether the true tree topology is useful
		* var-calls.vcf - gatk HaplotypeCaller calls
		* depth-and-het-filter.vcf - var-calls.vcf filtered to remove sites with too high DP and ExcessHet
		* repeat-filter.vcf - depth-and-het-filter.vcf filtered to remove repeatmasked regions
		* replicate-filter-only-variable.vcf - repeat-filter.vcf filtered to remove sites where any of the replicates of a sample disagree
		* replicate-filter-only-variable.fa - a FASTA alignment of the concatenated variable sites
		* replicate-filter-only-variable.nwk - a maximum likelihood tree of the FASTA alignment
 	* :open_file_folder: variant_table/ - table containing some extra data about each variant call
		* var_calls.tsv - table containing extra data about each variant call
		* chromosome_plot.pdf - plot of the first 11 chromosomes with variant locations indicated
 * :open_file_folder: results/
	* dng_callability.txt - analysis of simulated mutations to determine false negative rate
	* dng_variants.tsv - analysis of variants detected by the pipeline
	* filtered_bed_excluded.fa - FASTA alignment of concatenated variable sites determined by GATK
	* genomecov.txt - stats for calculating fraction of genome covered by a gene
	* MC_variants.bed - position of variants on the resistant branch
	* MC_closest_features.bed - gene annotations closest to the positions mutated on the resistant branch
	* MC_closest_genes.txt - genes in MC_closest_features.bed, only a subset of columns
	* MC_genes_phytozome_descriptions.tsv - phytozome descriptions for genes in MC_closest_genes.txt
 * :open_file_folder: scripts/
	* cigarpos.pl - helper script for induce_mutations.sh
	* clean_reads.sh - correct reads and remove reads with high-abundance kmers
	* create_consensus.sh - take a read alignment and a reference and produce a fasta consensus sequence
	* diploidify.py - take a multiple sequence alignment with IUPAC ambiguity codes representing heterozygous genotypes, do some checks, and turn it into one where each site represents one base
	* fig-em-compare.R - plots the true tree back to back with the inferred tree
	* fig-em-hist.R - plots the context histogram
	* filt_with_replicates.pl - filter a vcf using sample replication
	* gatkcaller.sh - an implementation of the GATK best practices workflow
	* induce_mutations.sh - generate sed scripts to modify reads from a simulated mutation
	* label_permutation.R - generate a random tree with maximal RF distance from the input tree
	* ngm_aligner.sh - align reads to a reference using NGM
	* optimize_parameters.R - find maximum likelihood DeNovoGear parameters
	* phytozome_query.py - query phytozome for gene descriptions
	* plot_chromosomes.R - plot the locations of variants on chromosomes
	* plot_tree.R - plot a newick tree
	* setup_mutations.sh - take a file of positions and genotypes and generate mutations with a TS:TV ratio of 2:1
	* slice-paired-reads-by-coverage.py - filter paired reads by a given kmer abundance cutoff
	* spiked_mutation_results.py - generates a table determining whether a mutation was detected given simulated mutations and variant calls
	* split_mutations_by_tree.sh - splits a table of proposed mutations to each sample such that a constant number of mutations are generated on each branch
	* variant_analysis.R - summarize data on the variants detected by the pipeline
	* vcf2fa.sh - helper script that takes a vcf and converts it to a diploidified FASTA alignment
	* vcf2tree.sh - helper script that takes a vcf and converts it to a newick tree view vcf2fa and RAxML
	* vcf_tab_to_fasta_alignment.pl - helper script to convert a vcf to a fasta alignment

Extras:
  * :open_file_folder: deprecated/
  * :open_file_folder: variant_analyses/

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
Variant calling can be performed using your favorite software.
We provide a script that works with short Illumina reads and uses GATK HaplotypeCaller to call variants.

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

Repeat Filters
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

Mutation Calling Given A .ped File with DeNovoGear
--------------------------------------------------
[DeNovoGear](https://github.com/denovogear/denovogear) is a suite of tools for detecting *de novo* mutations.
We can use it to find *de novo* mutations in our samples by encoding the sample structure as a pedigree
in a `.ped` file. In our case, the file looks like this:

```
##PEDNG v1.0
M       .       .       0       (((((M1a:0,M1b:0,M1c:0)M1:698,(M2a:0,M2b:0,M2c:0)M2:638)MD:75,(M3a:0,M3b:0,M3c:0)M3:796)MC:773,(M4a:0,M4b:0,M4c:0)M4:844)MB:372,((((M7a:0,M7b:0,M7c:0)M7:978,(M6a:0,M6b:0,M6c:0)M6:928)MG:111,(M8a:0,M8b:0,M8c:0)M8:1029)MF:112,(M5a:0,M5b:0,M5c:0)M5:1297)ME:358)MA:0;
```

To summarize, each line represents an individual. The 1st column is the individual ID, the 2nd and 3rd column represent the individual's father and mother.
The fourth column represents the individual's sex, where 0 indicates autosomal. The last column is a Newick tree listing each sample and a scaling factor
based on the physical distance between the samples. For more information on the PEDNG format, see [here](https://github.com/denovogear/denovogear#pedigree-file-format).

DeNovoGear can call variants from an alignment or pileup. It can also re-call sites provided in a VCF file.
We use DeNovoGear to call variants from HaplotypeCaller output with the `--standard-min-confidence-threshold-for-calling=0` parameter
set to benefit form HaplotypeCaller's haplotype reconstruction and DeNovoGear's variant calling model.

An example call using DeNovoGear this way is:

```bash
dng call --input=gatk.vcf --output=dng.vcf --ped=sampleM.ped --mu-somatic=1e-7 --theta=0.025
```

These calls can then be filtered using your favorite filtering method.

Estimating False Positive and False Negative Rates
--------------------------------------------------
False positive rates are difficult to estimate, but we can get a crude estimate by randomizing the assumed phylogeny,
rerunning the pipeline, and counting the number of new sites that are called. This estimates the false positive rate
by estimating how many noisy sites appear to be true variants by adding the tree to the model. Our pipeline produced
a mean of 0, but your mileage may vary.

We provide an R script `label_permutation.R` which shuffles the tree labels such that the new tree is maximally distant
from the given tree. The script takes no arguments but should be easily modifiable for use with other tree-like samples.
The following command will produce a pedigree file like the one above but with suffled labels, call DeNovoGear once more,
and remove the sites that were called with the original tree:

```bash
echo "##PEDNG v1.0" > shuffled.ped && echo -e "M\t.\t.\t0" | paste - <( label_permutation.R ) >> shuffled.ped #create a shuffled pedigree file
dng call --input=gatk.vcf --output=shuffled.vcf --ped=shuffled.ped --mu-somatic=1e-7 --theta=0.025 #call with shuffled tree
# ...
#apply filters here
# ...
bcftools view -T ^shuffled.vcf -Oz -o shuffled_no_original_calls.vcf.gz
```

The last step can be omitted, but you run the risk of including truly variable sites in your assessment.
Whether to include them or not depends on whether you would rather risk your rate being inflated or deflated.

The method we use to estimate false negative rates is to simulate mutations, re-run our pipeline, and
calculate the proportion of inserted mutations we ultimately recover.

TODO: talk about estimating FNR

Experimental Replication
========================

Makefiles
---------
We provide Makefiles in the analysis/ and data/ directories to allow for replication of our experiment and as an example of a pipeline.
We provide default variables, but they can be overwritten by setting the variable while calling `make`; for example, `make REF=myreference.fa`.

```bash
cd analysis
make
```

TODO:provide instruction on where you could swap in your reads to execute the same pipeline.
TODO:provide more examples

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

