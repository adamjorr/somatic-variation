#!usr/bin/perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use List::Util qw(any);
use Vcf;

my $vcffile ='/home/adam/eucalyptus/gatk/strict/filtered.vcf';
my $samfile ='/home/adam/eucalyptus/snap/snap.sam';

my $vcf = Vcf->new(file => $vcffile);
$vcf -> parse_header();
my %chrs;
while (my $line = $vcf -> next_line()){
	(my $id, my $loc, my $ref, my $alt) = (split("\t",$line))[0,1,3,4];
	my $snp = $ref . '/' . $alt;
	defined $chrs{$id} or $chrs{$id} = [];
	push(@{$chrs{$id}}, [$loc, $snp]);
}


open(my $samfh, "<", $samfile) or die $!;
print join("\t",qw(snp chromosome discosnploc gatkloc gt)) . "\n";
while(<$samfh>){
	chomp;
	next if substr($_,0,1) eq '@';
	(my $qname, my $rname, my $pos) = (split('\t'))[0,2,3];

	my $lunitiglen = (split('_',(split('\|',$qname))[6]))[-1];
	my $discogt = (split('_',(split('\|',$qname))[1]))[-1];
	my $bubpos = (split(':',(split('_',(split('\|',$qname))[1]))[-2]))[-1];
	my $discosnploc = $pos + $lunitiglen + $bubpos;

	next unless exists $chrs{$rname};
	my @gatklocs = @{$chrs{$rname}};
	#print(join('\t',($qname, $rname,$pos $snploc)));

	if (my @possible_loc_gts = grep {abs(${$_}[0] - $discosnploc) < 100} @gatklocs){
		for my $loc_gts (@possible_loc_gts){
			my @loc_gt = @{$loc_gts};
			if ($loc_gt[1] eq $discogt or reverse($loc_gt[1]) eq $discogt){
				print join("\t",((split('\|',$qname))[0], $rname, $discosnploc, $loc_gt[0], $loc_gt[1])), "\n";
			}
		}
	}

}








