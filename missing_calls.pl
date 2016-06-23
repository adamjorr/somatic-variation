#!usr/bin/env perl
#reads a VCF file from STDIN and prints the fraction of GTs containing a missing call.

use Vcf;
use List::MoreUtils qw(any);

my $vcf = Vcf->new(fh=>\*STDIN);
my $missing = 0;
my $total = 0;


$vcf->parse_header();
my @samples = $vcf->get_samples();



while(my $record = $vcf -> next_data_array()){
	my $fmt = $vcf->get_column($record,'FORMAT');
	my $gt_idx = $vcf->get_tag_index($fmt,'GT',':');

	if ($gt_idx == -1){
		warn "GT field not found for $record\nSkipping.";
		continue;
	}

	for my $sample (@samples){
		$total++;
		my $sam_column = $vcf->get_column($record,$sample);
		my $gt = $vcf->get_field($sam_column,$gt_idx);
		my @gts = $vcf->split_gt($gt);
		if (any {$_ eq '.'} @gts){
			$missing++;
		}
	}
}

print join("\t",qw(TOTAL MISSING PERCENT_MISSING)),"\n";
print join("\t",$total,$missing,$missing/$total),"\n";

exit;