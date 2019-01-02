# data/mutated

The following files come from the analysis/ make file
 * `e_mel_3_randsites.bed`
 * `e_mel_3_randsites.regions`
 * `overlapping_randsites.ad`

The following file contains the first 11 scaffold names, and currently lives in `/storage/datasets/Eucalyptus_melliodora/sampleM.regions`
 * sampleM.regions

Then produce a VCF to be mutated:
 * `dng call --lib-bias=2.473 --lib-error=0.002 --lib-overdisp=0.184 --mu=1.0e-08 --mu-somatic=0.0 --mu-library=0.0 --model='autosomal' --nuc-freqs='0.3,0.2,0.2,0.3' --ref-weight=4.93 --theta=0.247 --fasta='' --header='' --ped='sampleM_star.ped' --region='@sampleM.regions' --rgtag='LB' --sam-files='' --output='output.bcf' --min-qlen=0 --min-basequal=13 --min-mapqual=0 --normalize-somatic-trees=true --min-prob=0.0 overlapping_randsites.ad`

To make the mutated fastqs:
```bash
cd ../../analysis/false_negative_rate/
make mutatedreads
```
