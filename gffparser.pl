#!usr/bin/env perl

use strict;
use warnings;
use Bio::Tools::GFF;
use Bio::Tools::CodonTable;
use Getopt::Long;
use Number::Closest;
use List::MoreUtils qw/firstidx/;
use Vcf;

my $gff_file;
my $annotation_file;
my $flanks_file;
my $vcf_file;
my %locations;
my %names;
my $help = 0;
my $counter = 0;
my $grouped = 1;

GetOptions(	'gff|g=s' => \$gff_file,
			'annotation|a=s' => \$annotation_file,
			'flanks|f=s' => \$flanks_file,
			'vcf|v=s' => \$vcf_file,
			'replicates|r=i' => \$grouped,
			'help|h|?' => \$help);

pod2usage(1) if $help;
die "Cannot read $gff_file" unless -r $gff_file;
die "Cannot read $annotation_file" unless -r $annotation_file;
die "Cannot read $flanks_file" unless -r $flanks_file;
die "Cannot read $vcf_file" unless -r $vcf_file;

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
my %cds_locations;
while (my $feature = $gffio->next_feature() ){
	my $feat = $feature-> primary_tag();
	if ($feat eq 'gene' or $feat eq 'CDS'){
		my $chr = $feature -> seq_id();
		if ($feat eq 'gene'){
			$locations{$chr} = {} unless exists $locations{$chr};
			my $name = ($feature -> get_tag_values('Name'))[0];
			$locations{$chr} -> {$feature->start} = $name;
			$locations{$chr} -> {$feature->end} = $name;
		    # print join("\t", $feature-> seq_id(), $feature->start, $feature->end, $name, $annot{$name}->{'name'}) , "\n";
		    # $counter++;
		    # exit if $counter > 10;
		}
		else{
			my $phase = $feature->frame;
			#print $gffio->gff_string($feature),"\n";
			#print $feature->frame,"\n";
			$cds_locations{$chr} = {} unless exists $cds_locations{$chr};
			$cds_locations{$chr} -> {$feature->start} = $phase;
			$cds_locations{$chr} -> {$feature->end} = $phase;
		}
	}
}

my %locs;
my %cds_locs;
for my $chr (sort keys %locations){
	my @all_locs = sort {$a <=> $b} keys $locations{$chr};
	my @all_cds_locs = sort{$a <=> $b} keys $cds_locations{$chr};
	$cds_locs{$chr} = \@all_cds_locs;
	$locs{$chr} = \@all_locs;
}

my %flanks;
open (my $ffh, "<", $flanks_file) or die $!;
while (<$ffh>){
	chomp;
	(my $id, my $loc, my $flank) = split("\t");
	$flanks{$loc} = $flank;
}
close $ffh or die $!;


my $vcf = Vcf->new(file=>"$vcf_file");
my $vcfheader = $vcf->parse_header();
my @samples = $vcf -> get_samples();
my @sampleindex;
for (my $i = 0; $i < scalar(@samples); $i+=$grouped){
	push(@sampleindex, $vcf->get_column_index($samples[$i]));
}



print join("\t", "snp", 1 .. scalar(@sampleindex), "loc", "gene", "distance", "go_terms", "coding_mut", "substitution") , "\n";

while (my $x = $vcf->next_data_array()){
	my $chr = $$x[0];
	my $location = $$x[1];
	my $coordinate = $chr . ':' . $location;
	my $flank = $flanks{$coordinate};
	next unless $flank;
	#print join("\t",@{$locs{$chr}}), "\n";
	my @cds_locs = @{$cds_locs{$chr}};
	my $closefinder = Number::Closest->new(number => $location, numbers => $locs{$chr});
	my $closest = $closefinder -> find;
	my $closestgene = $locations{$chr}->{$closest};
	my $distance = abs($closest - $location);
	my $genename = $annot{$closestgene}->{'name'} || $closestgene;
	my $gocat = $annot{$closestgene}->{'go'};
	my $in_CDS;
	my $phase;
	my $old_aa;
	my $new_aa;
	my @muts = map{(split(':',$$x[$_]))[0]} @sampleindex;

	my $firstindex = List::MoreUtils::firstidx { $_ > $location } @cds_locs;
	if ($firstindex % 2 == 0){ #if 1st loc after snp is even
		#it's the start of a CDS and so is out of a CDS
		$in_CDS = 0;
		print join("\t", $flank, @muts, $coordinate, $genename, $distance, $gocat, $in_CDS,""),"\n";
	}
	else{#it's in a CDS
		$in_CDS = 1;
		$phase = $cds_locations{$chr} -> {$cds_locs[$firstindex]};
		my $firstbaseloc = $cds_locs[$firstindex - 1] + $phase;
		my $snp_codon_pos = (($location - $firstbaseloc) % 3);
		$flank =~ m/[ATCG]/;
		my $snp_flank_pos = $-[0];
		my $old_codon = substr($flank, $snp_flank_pos - $snp_codon_pos, 3);
		my $new_codon = $old_codon;
		substr($new_codon,$snp_codon_pos,1,$$x[4]);
		#print $new_codon, "\t", $old_codon, "\n";
		my $codontable = Bio::Tools::CodonTable->new();
		$old_aa = $codontable->translate($old_codon);
		$new_aa = $codontable->translate($new_codon);
		print join("\t", $flank, @muts, $coordinate, $genename, $distance, $gocat, $in_CDS, $old_aa . "/" . $new_aa),"\n";
	}
	#col 8 in gff. phase 0 means the first base is the first base of a codon, 1 means the 2nd base is the first base of a codon, 2 means the 3rd base is the first base of a codon.
}


$gffio -> close();
exit;

__END__

1. Sequence of flanking region around variant	x
2. 8 columns, each with a GT of a sample
3. Location on grandis	x
4. Nearest gene in grandis	x
5. Distance to nearest gene in grandis	x
6. GO categories	x
7. Mutation in CDS or not	x
8. AA change if in CDS

TODO:Some SNPs aren't getting their closest gene name and GO terms!!!