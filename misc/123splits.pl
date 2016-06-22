#!usr/bin/env perl
#123splits.pl
#A quick check of the 123 split.

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use List::MoreUtils qw(any all);
use Pod::Usage;
use Vcf;

my $seqio = Bio::SeqIO->newFh(-format => 'Fasta', -fh => \*STDIN);
my $seqout = Bio::SeqIO->newFh(-format => 'Fasta', -fh => \*STDOUT);
while(<$seqio>){
	my %sam_gt;
	my $id = $_ -> id();
	my @gts = grep {m/G[0-9]*.+/} split('\|',$id);
	@gts = map {(split(':',$_))[0]} @gts;

	%sam_gt = map{(split('_',$_))[0] => (split('_',$_))[1]} @gts;
	my @desired = qw(G1 G4 G7);
	my @others = qw(G10 G13 G16 G19 G22);
	if ($sam_gt{'G1'} eq $sam_gt{'G4'} and $sam_gt{'G4'} eq $sam_gt{'G7'}){
		if(all {$sam_gt{$_} ne $sam_gt{$desired[0]}} @others){
			print $seqout $_;
		}
	}
}

exit;