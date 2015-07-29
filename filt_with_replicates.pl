#!usr/bin/env perl
#perl filt_with_replicates.pl [-s | --strict] [-g | --grouped int] <vcfile.vcf >filtered.vcf

#Filters a VCF file by removing sites for a sample if all its replicates don't match.
#Using -s will remove the site for all samples if one of its samples' replicates don't match.

=head1 NAME

filt_with_replicates - Filter a VCF using replicate information.

=head1 SYNOPSIS

perl filt_with_replicates.pl [-s | --strict] [-g | --grouped int] [-h | -? | --help] <vcfile.vcf >filtered.vcf

Options:

=over 8

=item --help		Prints help information

=item --grouped		VCF samples aren't named, but in groups of int.

=item --strict		Use strict filtering

=back

=head1 OPTIONS

=over 8

=item B<--help>

Prints help information.

=item B<--grouped INT>

Use this option when the samples in the VCF aren't properly named OR when they are in a group.
This will assume that every INT genotypes in the VCF are replicates of one sample.
ONLY use this option when the replicates are in the correct order in the VCF.
For example, replicates in a VCF with a grouping of 3 might be: 1, 1, 1, 2, 2, 2, 3, 3, 3.

=item B<--strict>

This will strictly filter, removing the site completely if there is one replicate that doesn't match.
Normally, only the sample with the disagreement is removed, rather than all samples for that site.

=back

=head1 DESCRIPTION

Filters a VCF file by removing sites for a sample if all its replicates don't match.

=cut


use strict;
use warnings;
use Vcf;
use Getopt::Long;
use Pod::Usage;

my $strict = '';
my $grouped = 0;
my $help = 0;
GetOptions(	'strict|s' => \$strict,
			'grouped|g=i' => \$grouped,
			'help|h|?' => \$help);
pod2usage(1) if $help;


my $vcf = Vcf->new(fh=>\*STDIN);
$vcf->parse_header();
print $vcf -> format_header();

my @samples = $vcf->get_samples();
my %replicates;

if ($grouped){
	my $truesamples = scalar(@samples) / $grouped;
	my $samnum = 0;
	while (@samples){
		my @few;
		while(scalar(@few) < $grouped){
			push(@few,shift(@samples))
		}
		$replicates{$samnum++} = \@few;
	}
}

else{
	for my $sample (@samples){
		my @a;
		my $name = substr($sample,0,-1);
		$replicates{$name} ||= \@a;
		push(@{$replicates{$name}},$sample);
	}
}


if ($strict){
	strict_filter()
}

else{
	basic_filter()
}

exit;

#------------------------------

=head1 FUNCTIONS

=over

=item B<strict_filter()>

Applies strict filtering to the VCF. That is,
only outputs lines from the VCF where every
replicate of every sample has the same genotype
and there is a sample which differs from the others.

=cut

sub strict_filter{
	while (my $record = $vcf -> next_data_array()){
		my $stop = 0;
		my @allgts;
		for my $sample (keys %replicates){
			my @rep = @{$replicates{$sample}};
			my @genotype = map { (split(':',$vcf -> get_column($record, $_)))[0]} @rep;
			if (@genotype != grep{$_ eq $genotype[0]} @genotype){
				#not equal
				$stop = 1;
				last;
			}
			next if ($genotype[0] eq './.');
			push @allgts, @genotype;
		}
		next if $stop == 1;
		next if @allgts == grep{$_ eq $allgts[0]} @allgts; #Comment this line if you want ALL
		print join("\t",@{$record}) . "\n";
	}
}

=item B<basic_filter()>

Applies basic filtering to the VCF. That is,
only outputs lines from the VCF where every
replicate of some sample has the same genotype
and there is a sample which differs from the others.

=back

=cut

sub basic_filter{
	while (my $record = $vcf -> next_data_array()){
		my @allgts;
		for my $sample (keys %replicates){
			my @rep = @{$replicates{$sample}};
			my @genotype = map { (split(':',$vcf -> get_column($record, $_)))[0]} @rep;
			if (@genotype != grep{$_ eq $genotype[0]} @genotype){
				#not equal
				map { @{$record}[$vcf -> get_column_index( $_ )] = './.:0,0:0.00' } @rep;
				next;
			}
			next if ($genotype[0] eq './.');
			push @allgts, @genotype;
		}
		next if @allgts == grep{$_ eq $allgts[0]} @allgts; #Comment this line if you want ALL
		print join("\t",@{$record}) . "\n";
	}
}