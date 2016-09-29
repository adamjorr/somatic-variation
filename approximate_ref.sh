#!/usr/bin/env bash

#This script creates an approximation of a reference for
#an individual given sequencing reads from that
#individual and a reference to a closely-related species.
#This works by mapping the reads to the reference,
#creating a consensus sequence to use as a new reference,
#and repeating for a given number of iterations.

USAGE="Usage: $0 [-t THREADS] [-I ITERATIONS] [-d TMPDIR] -r reference.fa -o out.bam data/"

















