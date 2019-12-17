# Analyses of variants produced by the pipeline

## Overview
These analyses describe the positive control analyses in the paper.

These analysis scripts start with the output of the GATK variant detection pipeline in the `/results` folder of this repository, and through a series of R scripts use those results to assess whether the somatic mutations detected by that pipeline are a better match for the physical structure of the Eucalyptus melliodora individual than would be expected by chance.

## analysis0_iqtree_analyses.R
This script just removes any gapped sites from the original output from the GATK pipeline. Depending on how you ran this pipeline, its output might not have gaps in the first place, but it is sensible to run this script anyway just in case.

1. Make sure you have the R libraries `phangorn` and `phytools` installed. 

2. Get a binary of IQtree that works on your system and put it into the `positive_control_analysis` folder.

3. Change the two hard-coded file paths in the R script so that they match paths on your machine.

4. Use a terminal, CD to the `positive_control_analysis` directory, then run the R script from the terminal like this:

`Rscript analysis0_alignments.R`

The script will create a folder called `/analysis0_alignments`, then put in a copy of the original alignment, and then a version of that alignment with slightly edited names (e.g. M1a -> M1), and without gaps.



## analysis1_positive_control.R
This analysis asks simply: does the tree we infer from a stringently filtered dataset of somatic variants match the structure of the physical tree more closely than we would expect by chance?

1. Make sure you have teh R library `ggplot2` installed
2. Change the hardcoded paths at the top of the script so that they match the paths on your machine.
3. Run the R script like this:

`Rscript analysis1_positive_control.R`

Note that the analysis is likely to take ~10 minutes or more, because it is running 1000 bootstraps of a maxmimum likelihood analysis in IQ-tree. 

The script will produce:

	* output files from IQ-TREE
	* `best_trees.txv` which contains all of the maximum parsimony trees and the ML tree, along with the p-values of whether they are closer to the true tree than expected by chance
	* `positive_control_figure.pdf` which contains a visual representation of the data in `best_trees.tsv` 


## analysis2_compare_brlens.R
This analysis tests whether the physical and genetic branch lengths are correlated, using regression through the origin (the correct analysis, since our expectation is that a zero-length physical branch cannot have any mutations), and for comparison a regression not through the origin (simply to double check that the intercept is not signficantly different from zero, which would indicate a poor model fit for the expectation above). It is important to note that this analysis uses the physical tree structure for the inference.

1. Make sure you have teh R libraries `adephylo` and `reshape2` installed
2. Change the hardcoded paths at the top of the script so that they match the paths on your machine.
3. Run the R script like this:

`Rscript analysis2_compare_brlens.R`


The script will produce three files:

	* `brlen_correlation_figure_forced_through_origin` a ggplot2 plot of the data and the model forced through the origin
	* `brlen_correlation_figure__not_forced_through_origin` a ggplot2 plot of the data and the model not forced through the origin
	* `model_fits.txt` the fit of the linear models from R


