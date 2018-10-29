#!usr/bin/env perl
#perl vcfmasker.pl [-p 95] -g file.gff <file.vcf >masked.vcf
#Removes lines from a VCF that occur in a repeatmasked region with probability > -p.

use strict;
use warnings;
use Bio::Tools::GFF;
use Getopt::Long;
use Vcf;
use List::Util qw/any/;

my $pvalue = 95.0;
my $gff_file;
my $help;
GetOptions(	'gff|g=s' => \$gff_file,
			'probability|p=f' => \$pvalue,
			'help|h|?' => \$help);

pod2usage(1) if $help;
die "Cannot read $gff_file" unless -r $gff_file;

main();

sub main{
	my $maskr = load_gff($gff_file);
	my $vcfr = print_vcf_header();
	print_vcf($vcfr,$maskr);
	exit;
}

sub load_gff{
	my $gff_file = shift;
	my %mask;
	my $gffio = Bio::Tools::GFF->new(-file => $gff_file, -gff_version => 3);
	while (my $feature = $gffio->next_feature() ){
		$mask{$feature->seq_id()} ||= [];
		push(@{$mask{$feature->seq_id()}},$feature);
	}
	return \%mask
}

sub print_vcf_header{
	my $vcf = Vcf->new(fh=>\*STDIN);
	$vcf->parse_header();
	print $vcf -> format_header();
	return \$vcf;
}

sub print_vcf{
	my $vcfr = shift;
	my $maskr = shift;
	my $vcf = ${$vcfr};
	my %mask = %{$maskr};
	while (my $x = $vcf->next_data_array()){
		my $chr = $$x[0];
		my $loc = $$x[1];
		next unless defined($mask{$chr});
		print join("\t",@{$x}) . "\n" unless (any {$_ -> contains($loc)} @{$mask{$chr}});
	}
}

__END__

Todo: test

add documentation


Due to improperly formatted GFF file . . .

------------- EXCEPTION: Bio::Root::BadParameter -------------
MSG: ' 3.8' is not a valid score
VALUE:  3.8
STACK: Error::throw
STACK: Bio::Root::Root::throw /usr/share/perl5/Bio/Root/Root.pm:486
STACK: Bio::SeqFeature::Generic::score /usr/share/perl5/Bio/SeqFeature/Generic.pm:468
STACK: Bio::Tools::GFF::_from_gff3_string /usr/share/perl5/Bio/Tools/GFF.pm:628
STACK: Bio::Tools::GFF::from_gff_string /usr/share/perl5/Bio/Tools/GFF.pm:434
STACK: Bio::Tools::GFF::next_feature /usr/share/perl5/Bio/Tools/GFF.pm:395
...
--------------------------------------------------------------

To fix, use s/\t /\t/g 