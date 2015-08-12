#!usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Pod::Usage;
use Getopt::Long;

my $infile


my %locmap;
open(my $infh, "<", $infile) or die $!;
while (<$infh>){
	chomp;
	(my $id, my $matchid, my $startloc) = (split('\t',$_))[0,1,8]
	my $lunitiglen = split('_',(split('\|',$id))[4])[-1]
	$locmap{$id} = ($matchid, $startloc + $lunitiglen);
}

#Now needs VCF parsing stuff to compare with


my $seqio = Bio::SeqIO->newFh(-format => 'Fasta',)


