#!usr/bin/env perl

use strict;
use warnings;
use Bio::Tools::GFF;
use Getopt::Long;
use Number::Closest;
use List::MoreUtils;

my $gff_file;
my $annotation_file;
my %locations;
my %names;
my $help = 0;
my $counter = 0;

GetOptions(	'gff|g=s' => \$gff_file,
			'annotation|a=s' => \$annotation_file,
			'help|h|?' => \$help);

pod2usage(1) if $help;
die "Cannot read $gff_file" unless -r $gff_file;
die "Cannot read $annotation_file" unless -r $annotation_file;

my %annot;
open(my $annotfh, "<", $annotation_file) or die "$!";
while(<$annotfh>){
	(my $jgi_name, my $go_terms, my $gene_name) = (split("\t"))[1,9,11];
	$go_terms ||= '';
	$gene_name ||= '';
	$annot{$jgi_name} = {'go' => $go_terms, 'name' => $gene_name};
	#print join("\t",$jgi_name, $annot{$jgi_name}), "\n";
	#$counter++;
	#exit if $counter > 10;
}

# print $annot{'Eucgr.A00001'}, "\n";
# exit;

my $gffio = Bio::Tools::GFF->new(-file => $gff_file, -gff_version => 3);

while(my $feature = $gffio->next_feature()) {
	next unless $feature -> primary_tag() eq 'gene';
	my $chr = $feature -> seq_id();
	$locations{$chr} = {} unless exists $locations{$chr};

	my $name = ($feature -> get_tag_values('Name'))[0];
	$locations{$chr} -> {$feature->start} = $name;
	$locations{$chr} -> {$feature->end} = $name;
    # print join("\t", $feature-> seq_id(), $feature->start, $feature->end, $name, $annot{$name}->{'name'}) , "\n";
    # $counter++;
    # exit if $counter > 10;
}

for my $snp in @snps{
	my $chr;
	my $location;
	my @all_locs = sort keys $locations{$chr};
	my $closefinder = Number::Closest->new(number => $location, numbers => \@all_locs);
	my $closest = $closefinder -> find;
	my $closestgene = $locations{$chr}->{$closest};
	my $distance = abs($closest - $location);
	my $genename = $annot{$closestgene}->{'name'};
	my $gocat = $annot{$closestgene}->{'go'};
	my $ingene;
	if (first_index {$_ > $location} @all_locs % 2 == 0){ #if 1st loc after snp is even
		#it's the start of a gene and so is out of a gene
		$ingene = 0;
	}
	else{#it's in a gene
		$ingene = 1;
	}
}


$gffio -> close();
