#!usr/bin/env perl
#perl get_flanks.pl -v in.vcf <in.fasta >out.fasta

=head1 NAME

get_flanks - Get flanking sequences for a SNP from DiscoSNP++

=head1 SYNOPSIS

perl get_flanks.pl [-d | --discosnp snps.vcf] [-v | --vcf snps.vcf] [-r | --ref reference.fasta] [-h | -? | --help] <discosnp.fasta >out.fasta

Use perldoc get_flanks.pl for even more help.

Options:

=over 8

=item --help

Prints help information

=item --vcf FILE

Name of VCF file

=back

=head1 OPTIONS

=over 8

=item B<--help>

Prints help information.

=item B<--vcf FILE>

This is the VCF which will be read to determine which SNPs to extract from the
input fasta file. This should be the file output by the DiscoSNP++ run_VCF_creator tool,
or some filtered subset of that VCF. B<get_flanks> is designed for SNPs and not INDELS or
other types of DiscoSNP calls, so your mileage may vary when attempting to use a VCF
with variants other than SNPs.

=back

=head1 DESCRIPTION

Reads DiscoSNP++-generated VCF files or some subset thereof and the accompanying
fasta file (generated with the -T option) and prints a fasta file with only the SNPs
specified in the VCF. The output fasta will have the SNP capitalized with all surrounding
sequences in lower case. Note that this is different from the default DiscoSNP++ output,
which prints the entire bubble in upper case and the surrounding contigs in lowercase.
Also note that only the first SNP in a bubble will be considered.

=cut

use strict;
use warnings;
use Vcf;
use Bio::SeqIO;
use Getopt::Long;
use List::MoreUtils qw(any);
use Pod::Usage;

my $vcffile;
my $help = 0;
my $reffile;
my $discofile;
GetOptions(	'discosnp|d=s' => \$discofile,
			'vcf|v=s' => \$vcffile,
			'ref|r=s' => \$reffile,
			'help|h|?' => \$help);
pod2usage(1) if $help;
die "Cannot read $discofile" unless -r $discofile;

my %chr_seq = chr_ids($reffile);
my @dids = dids_from_vcf($discofile);
print_fasta(@dids);
exit;


# -----------------------

sub chr_seq{
	my $reffile = shift;
	my %chrids;
	my $seqio = Bio::SeqIO->newFh(-format => 'Fasta', -file => $reffile);
	while(<$seqio>){
		$chrids{$_ -> id()} = $_ -> seq();
	}
	return %chrids;
}

sub dids_from_vcf{
	my $vcffile = shift;
	my %chr_seq = shift;
	my %ids;
	my $vcf = Vcf->new(file=>"$vcffile");
	$vcf->parse_header();
	while (my $x = $vcf -> next_data_array()){
		my $chr = $$x[0];
		if (any {$_ == $chr} keys %chr_seq){
			$ids{$$x[2]} = {'chr' => $chr, 'pos' => $$x[1]};
		}
		else{
			$ids{$$x[2]} = {'chr' => '', 'pos' => ''};
		}
	}
	return %ids;
}

sub flanking_ref{
	my $chr = shift;
	my $pos = shift;
	my $flanksize = shift;
	my %chrids = shift;

	my $seq = $chrids{$chr};
	my $startloc;
	if ($flanksize > $pos + 1){
		$startloc = 0;
	}
	else{
		$startloc = $pos + 1 - $flanksize;
	}

	return lc(substr($seq, $startloc, $flanksize)) . uc(substr($seq,$pos,1)) . lc(substr($seq, $pos + 1, $flanksize));

}

sub flanking_disco{
	my %did_seq = shift;
	my $id = shift;
	my $pos = shift;

	my $seq = $did_seq{$id};
	return lc(substr($seq,0,$pos)) . uc(substr($seq,$pos,1)) . lc(substr($seq,$pos + 1));
}

sub discoid_seq{
	my %did_seq;
	my $seqio = Bio::SeqIO->newFh(-format => 'Fasta', -fh => \*STDIN);
	while (<$seqio>){
		my $full_id = $_ -> id();
		my $id = (split('_',(split('\|',$full_id))[0]))[-1];
		$did_seq{$id} = $_ -> seq();
	}
}

sub othervcflocs{
	my $vcffile = shift;
	my %chr_seq = shift;
	my %ids;
	my $counter = 0;

	my $vcf = Vcf->new(file=>"$vcffile");
	while (my $x = $vcf -> next_data_array()){
		my $chr = $$x[0];
		my $pos = $$x[1];
		$ids{$counter} = {'chr' => $chr, 'pos' => $pos}
	}
}






sub print_fasta{
	my @ids = @_;
	my $seqio = Bio::SeqIO->newFh(-format => 'Fasta', -fh => \*STDIN);
	my $seqout = Bio::SeqIO->newFh(-format => 'Fasta', -fh => \*STDOUT);
	while (<$seqio>){
		my $full_id = $_ -> id();
		(my $name, my $full_snp) = (split('\|',$full_id))[0,1];
		if (any {$_ eq $name} @ids){
			my $snp = (split(',',$full_snp))[0];
			my @split = split('_',$snp);
			(my $num, my $pos) = split(':',$split[1]);
			(my $ref, my $alt) = split('\/',$split[2]);
			my $seqstr = $_ -> seq();
			my ($uppers) = $seqstr =~ m/[A-Z]+/g;
			my $upperstart = index($seqstr,$uppers);
			$seqstr = lc($seqstr);
			substr($seqstr,$upperstart + $pos,1) = uc(substr($seqstr,$upperstart + $pos,1));
			$_ -> seq($seqstr);
			print $seqout $_;
		}
	}
}

__END__

Output 3 files: unique DISCOSNP variants, unique GATK variants, and both variants.
DISCOSNP variants should have a DISCOSNP ID, LOCATION, FLANKS.
GATK variants should have a LOCATION, FLANKS.





