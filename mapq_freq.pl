#!usr/bin/env perl
#Returns a tab-delimited file of frequencies of a list of numbers given, one on each line.
#For example, samtools view FILE | cut -f5 | mapq_freq.pl > stats.txt

use strict;
use warnings;

my %freqs;
while(<>){
	chomp;
	$freqs{$_} ||= 0;
	$freqs{$_} += 1;
}

for my $k (sort keys %freqs){
	print $k . "\t" . $freqs{$k} . "\n";
}

exit;