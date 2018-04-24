# Analyses of variants produced by the pipeline

## Overview
These analysis scripts start with the output of the variant detection pipeline in the `/results` folder of this repository, and through a series of R scripts use those results to measure various quantities about somatic mutation.

## analysis0_iqtree_analyses.R
To run this script, you will need to have a binary of IQtree that works on your system and put it into the `analysis0_iqtree` folder.

This script writes out 5 alignments into the `variant_analyses/alignments` folder. Each alignment permits a different number of gaps in each column:

* `aln0.fa`: no gaps in any column
* `aln1.fa`: up to 1 gaps per column
* `aln2.fa`: up to 2 gaps per column
* `aln3.fa`: up to 3 gaps per column
* `aln4.fa`: up to 4 gaps per column
* `aln5.fa`: up to 5 gaps per column

It then runs IQtree on each alignment, with default settings and 250 bootstrap replicates.

## analysis1_positive_control.R
This analysis asks simply: does the tree we infer from a stringently filtered dataset of somatic variants match the true structure of the physical tree more closely than we would expect by chance?


