use Bio::Cigar;

# this script takes in two params:
# -c=<cigar string>
# -p=<unadjusted position of SNP in an alignment in the BAM>
# example:  perl -s cigarpos.pl -c=10S50M2D10M -p=30

# print "translating BAM alignment pos $p to FASTQ read pos\n";

my $cigar = Bio::Cigar->new($c);

#CORE::say "Query length is ", $cigar->query_length;
#CORE::say "Reference length is ", $cigar->reference_length;
my ($qpos, $op) = $cigar->rpos_to_qpos($p);
#CORE::say "query pos at ref pos ",$p," is ", $qpos;
#CORE::say "Alignment operation at ref pos ",$p," is ", $op;
if($qpos ne "") {
	print $qpos;
}
else {
	print -1;
}
