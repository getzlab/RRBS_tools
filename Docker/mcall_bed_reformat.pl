#!/opt/miniconda/bin/perl
  
use strict;
use warnings;



## quick&dirty: takes just one argument
unless (@ARGV==1) {
	print "v06.11.2018";
    exit;
}

my $bed = shift;

(my $bed_modified = $bed) =~ s/(.+)\.bed$/${1}.pre-bb.bed/;

open(my $fh_r, "<", $bed)
    or die "[ERROR] Cannot read BED file '$bed': $!\n";

open(my $fh_w, ">", $bed_modified)
    or die "[ERROR] Cannot write BED file '$bed_modified': $!\n";

## original
## chr10    3196577 3196579 0.6 5   3   B   G   +   2   0   -   3   3   GGCGC
## chr10    3463843 3463845 0   5   0   B   G   +   4   0   -   1   0   GTCGG
## chr10    7229647 7229649 0.273   11  3   B   G   +   10  2   -   1   1   TGCGT
##
## new format (pre-bb)
## chr10    3196577 3196579 '60%[5]'    600 +
## chr10    3463843 3463845 '0%[5]' 0   +
## chr10    7229647 7229649 '27%[11]'   273 +

my @rec;
my $str;

while (my $line = readline($fh_r)) {

    ## #chrom start end ratio totalC
    next if $line =~ /^#/;

    @rec = split(/\t/,$line);

    ## chr10 7230389 7230391 0.333 12 <..>
    next unless $rec[4]>4;

    $str  = "$rec[0]\t";
    $str .= "$rec[1]\t";
    $str .= "$rec[2]\t";
    $str .= "\'" . int($rec[3]*100) . "%[" . $rec[4] . "]\'" . "\t";
    $str .= int($rec[3]*100) . "\t";
    $str .= "\t+\n";

    print $fh_w $str;
}
close $fh_r;
close $fh_w;

print "$bed_modified";
