A Pipeline For the Detection of Somatic Variants
================================================
Herein we describe our pipeline for detecting somatic mutants.

Software Prerequisites
----------------------

Input File Requirements
-----------------------

Experimental Overview
---------------------
(Explain the conditions surrounding the experiment here)


See the section on replicating our experiment for more information.



Creating A Reference From a Similar Species with Iterative Mapping
==================================================================
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
Here will be a script once I make sure the method works
```





Experimental Replication
========================














