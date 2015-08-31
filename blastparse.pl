#!usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Pod::Usage;
use Getopt::Long;
use List::Util qw(any);
use Vcf;

my $blastfile = '/home/adam/eucalyptus/blast/results.tab';
my $vcffile ='/home/adam/eucalyptus/gatk/strict/filtered.vcf';

my $vcf = Vcf->new(file => $vcffile);
$vcf -> parse_header();
my %chrs;
while (my $line = $vcf -> next_line()){
	(my $id, my $loc, my $ref, my $alt) = (split("\t",$line))[0,1,3,4];
	my $snp = $ref . '/' . $alt;
	defined $chrs{$id} or $chrs{$id} = [];
	push(@{$chrs{$id}}, [$loc, $snp]);
}

open(my $infh, "<", $blastfile) or die $!;
print join("\t" , qw(disco_id chr disco gatkloc gt)) , "\n";
while (<$infh>){
	chomp;
	(my $id, my $matchid, my $qstart, my $alnlen, my $sstart) = (split('\t',$_))[0,1,4,5,7];
	my $lunitiglen = (split('_',(split('\|',$id))[6]))[-1];
	my $discogt = (split('_',(split('\|',$id))[1]))[-1];

	next if $qstart > $lunitiglen;
	next if $qstart + $alnlen < $lunitiglen;
	my $discosnploc = $sstart + $lunitiglen - $qstart;

	next unless exists $chrs{$matchid};
	my @gatklocs = @{$chrs{$matchid}};

	if (my @possible_loc_gts = grep {abs(${$_}[0] - $discosnploc) < 1000} @gatklocs){
		for my $loc_gts (@possible_loc_gts){
			my @loc_gt = @{$loc_gts};
			if ($loc_gt[1] eq $discogt or reverse($loc_gt[1]) eq $discogt){
				print join("\t",((split('\|',$id))[0], $matchid, $discosnploc, $loc_gt[0], $loc_gt[1])), "\n";
			}
		}
	}
}

##Make a check that it's the same SNP



