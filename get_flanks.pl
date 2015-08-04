#!usr/bin/env perl
#perl get_flanks.pl -v in.vcf <in.fasta >out.fasta

use strict;
use warnings;
use Vcf;
use Bio::SeqIO;
use Getopt::Long;


my $vcf = Vcf->new(fh=>\*STDIN);
$vcf->parse_header();

my $seqio = Bio::SeqIO->newFh(-format => 'Fasta', -fh => \*STDIN);
while (<$seqio>){
	my $name = (split('\|',$_ -> id()))[0];
	print $name , "\n"
}