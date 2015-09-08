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
my $flanksize = 200;
my $reffile;
my $discofile;
GetOptions(	'discosnp|d=s' => \$discofile,
			'vcf|v=s' => \$vcffile,
			'ref|r=s' => \$reffile,
			'help|h|?' => \$help);
pod2usage(1) if $help;
die "Cannot read $discofile" unless -r $discofile;

print_output();
exit;


# -----------------------

sub chr_seq{
	my $reffile = shift;
	my %chrids;
	my $seqio = Bio::SeqIO->newFh(-format => 'Fasta', -file => $reffile);
	while(<$seqio>){
		$chrids{$_ -> id()} = $_ -> seq();
	}
	return \%chrids;
}

sub dids_from_vcf{
	my $vcffile = shift;
	my $chr_seq_ref = shift;
	my %chr_seq = %$chr_seq_ref;
	my %ids;
	my $vcf = Vcf->new(file=>"$vcffile");
	$vcf->parse_header();
	while (my $x = $vcf -> next_data_array()){
		my $chr = $$x[0];
		if (any {$_ eq $chr} keys %chr_seq){
			$ids{$$x[2]} = {'chr' => $chr, 'pos' => $$x[1]};
		}
		else{
			$ids{$$x[2]} = {'chr' => '', 'pos' => $$x[1]};
		}
	}
	return \%ids;
}

sub flanking_ref{
	my $chr = shift;
	my $pos = shift;
	my $flanksize = shift;
	my $chridsref = shift;
	my %chrids = %$chridsref;

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
	my $did_seqref = shift;
	my %did_seq = %$did_seqref;
	my $id = shift;
	my $pos = shift;

	my $seq = $did_seq{$id};
	die "Can't find seq for $id!" unless $seq;
	return lc(substr($seq,0,$pos)) . uc(substr($seq,$pos,1)) . lc(substr($seq,$pos + 1));
}

sub discoid_seq{
	my %did_seq;
	my %did_loc;
	my $seqio = Bio::SeqIO->newFh(-format => 'Fasta', -fh => \*STDIN);
	while (<$seqio>){
		my $full_id = $_ -> id();
		my $id = (split('_',(split('\|',$full_id))[0]))[-1];
		$did_seq{$id} = $_ -> seq();

		my $loc = (split(':',(split('_',(split('\|',$full_id))[1]))[1]))[1] + (split('_',(split('\|',$full_id))[6]))[-1];
		$did_loc{$id} = $loc;
	}
	return (\%did_seq , \%did_loc);
}

sub othervcflocs{
	my $vcffile = shift;
	my $chr_seq_ref = shift;
	my %chr_seq = %$chr_seq_ref;
	my %ids;
	my $counter = 0;

	my $vcf = Vcf->new(file=>"$vcffile");
	$vcf -> parse_header();
	while (my $x = $vcf -> next_data_array()){
		my $chr = $$x[0];
		my $pos = $$x[1];
		$ids{$counter} = {'chr' => $chr, 'pos' => $pos};
		$counter++;
	}
	return \%ids;
}

sub find_mutual_snps{
	#(%left, %both, %right)
	my $did_locr = shift;
	my $otherid_locr = shift;
	my %did_loc = %$did_locr;
	my %otherid_loc = %$otherid_locr;
	my %left;
	my %both;
	my %right;

	for my $key (sort keys %did_loc){
		my $dloc = $did_loc{$key};
		if (any {$otherid_loc{$_}->{'chr'} eq $dloc->{'chr'} && $otherid_loc{$_}->{'pos'} == $dloc->{'pos'}} keys %otherid_loc){
			$both{$key} = $did_loc{$key};
		}
		else{
			$left{$key} = $did_loc{$key};
		}
	}

	for my $key (sort keys %otherid_loc){
		if (! any {%{$did_loc{$_}} == %{$otherid_loc{$key}}} ){
			$right{$key} = $otherid_loc{$key};
		}
	}
	return (\%left, \%both, \%right);
}

sub print_output{
	my $chr_seq_ref = chr_seq($reffile);
	my $didsr = dids_from_vcf($discofile, $chr_seq_ref);
	my $otherids = othervcflocs($vcffile, $chr_seq_ref);
	(my $donlyr, my $bothr, my $otheronlyr) = find_mutual_snps($didsr, $otherids);
	my %donly = %$donlyr;
	my %both = %$bothr;
	my %otheronly = %$otheronlyr;

	(my $dseqsr, my $dlocr) = discoid_seq();
	my %dseqs = %$dseqsr;
	my %dlocs = %$dlocr;
	my %discoflanks;
	my %otherflanks;
	open(my $discofh, ">", 'discoflanks.txt') or die $!;
	open(my $bothfh, ">", 'bothflanks.txt') or die $!;
	open(my $otherfh, ">", 'gatkflanks.txt') or die $!;
	print{$discofh}(join("\t",qw'DISCOSNP_ID LOC FLANKS'),"\n");
	print{$bothfh}(join("\t",qw'DISCOSNP_ID LOC FLANKS'),"\n");
	print{$otherfh}(join("\t",qw'LOC FLANKS'),"\n");


	for my $key (sort keys %donly){
		my $pos = $donly{$key}->{'pos'};
		my $flanks = flanking_disco($dseqsr, $key, $dlocs{$key});
		my $chr = $donly{$key}->{'chr'};
		if ($chr){
			print{$discofh}(join("\t",$key, $chr . $pos, $flanks),"\n");	
		}
		else{
			print{$discofh}(join("\t",$key, '', $flanks),"\n");	
		}
		
	}
	close $discofh;

	for my $key (sort keys %both){
		my $pos = $both{$key}->{'pos'};
		my $flanks = flanking_disco($dseqsr, $key, $dlocs{$key});
		print{$bothfh}(join("\t",$key, $both{$key}->{'chr'} . $pos, $flanks),"\n");
	}
	close $bothfh;

	for my $key (sort keys %otheronly){
		my $pos = $otheronly{$key}->{'pos'};
		my $chr = $otheronly{$key}->{'chr'};
		my $flanks = flanking_ref($chr, $pos, $flanksize, $chr_seq_ref );
		print{$otherfh}(join("\t",$chr . $pos, $flanks),"\n");
	}
	close $otherfh;

}

# ------------------
# DEPRECATED
# ------------------
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

Fix disco flanks; the position is not correctly given.
Pos is inaccurate if the SNP is mapped




