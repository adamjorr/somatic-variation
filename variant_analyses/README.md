# Analyses of variants produced by the pipeline

## Overview
These analysis scripts start with the output of the variant detection pipeline in the `/results` folder of this repository, and through a series of R scripts use those results to measure various quantities about somatic mutation.

## analysis0_write_alignments.R
This script writes out 7 alignments into the `variant_analyses/alignments` folder. Each alignment permits a different number of gaps (coded as Ns) in each column:

* `aln0.fa`: no gaps in any column
* `aln1.fa`: up to 1 gaps per column
* `aln2.fa`: up to 2 gaps per column
* `aln3.fa`: up to 3 gaps per column
* `aln4.fa`: up to 4 gaps per column
* `aln5.fa`: up to 5 gaps per column

With 5 gaps per column, that column will have just 3 characters (we have 8 branches in total), so has the potential to be phylogenetically informative. If we allow more gaps, there is no chance for a column to be phylogenetically informative, so we don't use it.

## analysis1_positive_control.R
This analysis asks simply: does the tree we infer from a stringently filtered dataset of somatic variants match the true structure of the physical tree more closely than we would expect by chance?


