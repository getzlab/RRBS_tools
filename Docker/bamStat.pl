#!/usr/local/bin/perl

## INFO : I dropped 'Bio::DB::Sam' approach in favour of 'samtools --threads INT' as this works 
##        slightly faster. And it is more familiar :-) - 2018-08-15

## $Header: /home/klages/bin/bamStat.pl,v 1.1 2014/10/20 13:28:01 klages Exp klages $
## $Log: bamStat.pl,v $
## Revision 1.1  2014/10/20 13:28:01  klages
## Initial revision

## H0VYCAGXX:1:8:4880407:0 99    chr2_part1  6757457  60 134M15S  = 6757457   133 AAAA[..]GCAC <AAA[..]FF.F SM:i:796 AS:i:1691   RG:Z:0   NM:i:0   BC:Z:none
## H0VYCAGXX:1:8:4880407:0 147   chr2_part1  6757457  60 15S134M  = 6757457  -133 ACGC[..]CCGG FFFF[..]AAAA SM:i:798 AS:i:1691   RG:Z:0   NM:i:0   BC:Z:none
## H0VYCAGXX:1:9:1230537:0 147   chr2_part1  6757477  60 119S30M  = 6757328  -178 AAAA[..]GGTG .FFF[..]AAAA SM:i:798 AS:i:1694   RG:Z:0   NM:i:0   BC:Z:none
## H0VYCAGXX:2:3:1366716:0 147   chr2_part1  6757498  60  62S87M  = 6757349  -235 TCAC[..]CCAG FFFF[..]AAAA SM:i:795 AS:i:1692   RG:Z:1   NM:i:0   BC:Z:none
## H0VYCAGXX:1:8:3961270:0 83    chr2_part1  6757514  60  1S148M  = 6757406  -255 AACG[..]GTGG .FFF[..]A.AA SM:i:798 AS:i:1694   RG:Z:0   NM:i:0   BC:Z:none
## H0VYCAGXX:1:3:4185308:0 99    chr2_part1  6757522  60    149M  = 6757671   252 GAAC[..]CACT <AAA[..]F<F< SM:i:781 AS:i:1675   RG:Z:0   NM:i:1   BC:Z:none
##
## NM:i : Edit distance to the reference, including ambiguous bases but excluding clipping (==M)
##
## M alignment match (can be a sequence match or mismatch)
## I insertion to the reference
## D deletion from the reference
## N skipped region from the reference
## S soft clipping (clipped sequences present in SEQ)
## H hard clipping (clipped sequences NOT present in SEQ)
## P padding (silent deletion from padded reference)
## = sequence match (E)
## X sequence mismatch
##
## The bitwise FLAG (position 2 in SAM):
## -------------------------------------
## Bit    Description
## 0x001  template having multiple segments in sequencing (PE in sequencing)
## 0x002  each segment properly aligned according to the aligner (mapped as proper pair)
## 0x004  segment unmapped
## 0x008  next segment in the template unmapped
## 0x010  SEQ being reverse complemented
## 0x020  SEQ of the next segment in the template being reversed
## 0x040  the first segment in the template
## 0x080  the last segment in the template
## 0x100  secondary alignment
## 0x200  not passing quality controls
## 0x400  PCR or optical duplicate
## 0x800  supplementary alignment
##
##
## For each read/contig in a SAM file, it is required that one and only one line
## associated with the read satisfies ‘FLAG & 0x900 == 0’. This line is called
## the primary line of the read.
##
## * Bit 0x100 marks the alignment not to be used in certain analyses when the tools
## in use are aware of this bit. It is typically used to flag alternative mappings
## when multiple mappings are presented in a SAM.
##
## * Bit 0x800 indicates that the corresponding alignment line is part of a chimeric
## alignment. A line flagged with 0x800 is called as a supplementary line.
##
## * Bit 0x4 is the only reliable place to tell whether the read is unmapped.
## If 0x4 is set, no assumptions can be made about RNAME, POS, CIGAR, MAPQ, bits 0x2,
## 0x10, 0x100 and 0x800, and the bit 0x20 of the previous read in the template.
##
## * If 0x40 and 0x80 are both set, the read is part of a linear template, but it is
## neither the first nor the last read. If both 0x40 and 0x80 are unset, the index of
## the read in the template is unknown. This may happen for a non-linear template or
## the index is lost in data processing.
##
## * If 0x1 is unset, no assumptions can be made about 0x2, 0x8, 0x20, 0x40 and 0x80.
##
## taken from:
## "Sequence Alignment/Map Format Specification"
## "The SAM/BAM Format Specification Working Group"
## "12 Sep 2014"
##

use strict;
use warnings;
use Getopt::Long;
use Statistics::Descriptive;

my %config;

## presets
$config{threads} = 4; # More than 4 does not give any benefit.
$config{verbose} = 0;
$config{chunk}   = 0; # how many seqs to analyze? 0=all
$config{alnType} = 'BAM';

GetOptions ( \%config,

   'sam=s',
   'bam=s',
   'verbose',
   'chunk=i',
   'threads=i',
   'version'

) or die "No! No. You don't understand. There was something in the woods, David... and I think it's in here with us... now.\n";

## simple version string
if($config{version}) {
	print "2018-08-15";
    exit;
}


## no buffering
$|++ if($config{verbose});


## (1) Read either BAM or SAM file.
my $fh;
if(exists $config{bam}) {
   open($fh,"samtools view -F 0x100 -F 0x800 --threads $config{threads} $config{bam} |")
      or die "Cannot read BAM file '$config{bam}': $!\n";
   $config{alnFile} = $config{bam};
}
elsif(exists $config{sam}) {
   open($fh,"<",$config{sam})
      or die "Cannot read SAM file '$config{sam}': $!\n";
   $config{alnType} = 'SAM';
   $config{alnFile} = $config{sam};
}
else {
   die "Please supply either BAM or SAM file.\n";
}


## these items are extracted drom SAM
my ($flag,$cigar,$sequence,$editDistance);

## pre-init cigar ops
my %cigar = map { $_ => 0 } qw(M I D N S H P X =);

my ($len,$type,$posInCigar); # split cigar
my $clip5p      = 0;  # clipped 5' bases
my $clip3p      = 0;  # clipped 3' bases
my $readCount   = 0;  # all reads with CIGAR available
my $readsUnmp   = 0;  # 0x004, unmapped
my $mappedPairs = 0;  # correctly mapped pairs
my $dupReads    = 0;  # duplicates
my $readLength;       # just in case reads were hardclipped
my $totalLength;      # total amount of bases of valid CIGAR sequences
my $opLength;         # sum of lengths of the M/I/S/=/X
my (@clip5p,@clip3p); # holds lengths of cut down bases, to calculate average length removed at either end
my @cigarGroups;

## timing
my $timeStart    = time();
my $timeEnd      = 0;

print STDERR "Reading $config{alnType} Alignment File '$config{alnFile}'.\n";

while($_=readline($fh)) {

   #next if(/^@/); # skip SAM fields

   ## iSAAC BAM lines; NM is not always written by other aligners. So as long as we don't use the edit distance for any statistics, we do not try to parse NM.
   ## H0VYCAGXX:1:8:3961270:0 83    chr2_part1  6757514  60  1S148M  = 6757406  -255 AACG[..]GTGG .FFF[..]A.AA SM:i:798 AS:i:1694   RG:Z:0   NM:i:0   BC:Z:none
   ## H0VYCAGXX:1:3:4185308:0 99    chr2_part1  6757522  60    149M  = 6757671   252 GAAC[..]CACT <AAA[..]F<F< SM:i:781 AS:i:1675   RG:Z:0   NM:i:1   BC:Z:none
   ##($flag,$cigar,$sequence,$editDistance)=$_=~/^[^\t]+\t([^\t]+)\t(?:[^\t]+\t){3}([^\t]+)\t(?:[^\t]+\t){3}([^\t]+)\t(?:[^\t]+\t){2,4}NM:i:(\d+)\t/;
   ($flag,$cigar,$sequence)=$_=~/^[^\t]+\t([^\t]+)\t(?:[^\t]+\t){3}([^\t]+)\t(?:[^\t]+\t){3}([^\t]+)/;

   ##                  Rate         split       pattern
   ## split()      940219/s           --         -43%
   ## pattern =~  1637003/s          74%           --

   ## skip unmapped reads
   if($flag & 0x004) {
      $readsUnmp++;
      next;
   }

   ## correctly mapped read pairs
   if(($flag & 0x001) && ($flag & 0x002) && ($flag & 0x040) && !($flag & 0x080)) {
      $mappedPairs++;
   }

   ## duplicated reads
   if($flag & 0x400) {
      $dupReads++;
   }

   $readCount++;
   $readLength   = length($sequence);
   $totalLength += $readLength;

   ## NOTE: def=0, supply on command line
   last if( ($readCount+$readsUnmp) == $config{chunk});

   if($config{verbose}) {
      print STDERR "." if(($readCount%1_000_000)==0);

      if(($readCount%10_000_000)==0) {
         print STDERR "[@{[commify($readCount)]} \@ @{[commify(calculateProcessRate($readCount,$timeStart))]} reads/s]";
      }
   }

   ## put spaces around OPs
   #$cigar=~s/(\d+[MIDNSHPX=])/ $1/g;

   ## remove very first whitespace
   #$cigar=~s/^\s+//;

   ## split CIGAR into OPs
   #@cigarGroups = split(' ', $cigar);

   @cigarGroups = split(/(\d+[MIDNSHPX=])/,$cigar);

   ## if read is mapped reversed, reverse CIGAR accordingly
   if($flag & 0x010) {
      @cigarGroups = reverse(@cigarGroups);
   }

   $posInCigar = 1; # 1=left end, >1 = not left end :-)
   $opLength   = 0;

   foreach my $op (@cigarGroups) {

      next unless ($op =~ /^\S+$/);

      ## M149 -> M, 149
      ($len,$type)=$op=~/(\d+)([MIDNSHPX=])/;

      ## sum up all itmes
      $cigar{$type} += $len;

      ## Sum of lengths of the M/I/S/=/X operations shall equal the length of $sequence
      if($type =~ /[MISX=]/) {
         $opLength += $len;
      }

      ## end clipping
      ## There are some cases like:
      ## a) 1M148S -> 3' clipping
      ## b) 148S1M -> 5' clipping
      if($type eq 'S' && $posInCigar==1) {
         $clip5p += $len;
         push(@clip5p,$len);
      }
      elsif($type eq 'S') {
         $clip3p += $len;
         push(@clip3p,$len);
      }
      $posInCigar++;
   }
   print STDERR "*** $cigar does not match sequence length. CIGAR=$opLength READ=$readLength ***\n" unless($opLength==$readLength);
}

## calculate final/average processing rate before statistics ops
my $processingRate = calculateProcessRate($readCount,$timeStart);

print STDERR "[@{[commify($readCount)]}]\n" if($config{verbose});

## Notbremse
unless($readCount) {
   print "No mapped reads found. Exit.\n";
   exit;
}

## Calculate some fractions
my %frac;
foreach my $op (keys %cigar) {
   $frac{$op} = sprintf("%.1f",$cigar{$op}*100/$totalLength);
}

$frac{validReads}    = sprintf("%.1f",$readCount*100/($readCount+$readsUnmp));
$frac{unmappedReads} = sprintf("%.1f",$readsUnmp*100/($readCount+$readsUnmp));
$frac{mappedPairs}   = sprintf("%.1f",$mappedPairs*100/$readCount);
$frac{dupReads}      = sprintf("%.1f",$dupReads*100/$readCount);
$frac{clipped5p}     = sprintf("%.1f",$clip5p*100/$totalLength);
$frac{clipped3p}     = sprintf("%.1f",$clip3p*100/$totalLength);

## Some stats for 5' end of read
my $stat5p = Statistics::Descriptive::Full->new();
$stat5p->add_data(@clip5p);

## The same for 3' end
my $stat3p = Statistics::Descriptive::Full->new();
$stat3p->add_data(@clip3p);

my $strLen = length($totalLength) + 3;
my $modalLength;

my $summary = "\n";
$summary .= "Input Statistics\n";
$summary .= "----------------\n";
$summary .= "\n";
$summary .= "  BAM File ......................: " . $config{bam} . "\n";
$summary .= "  Number of Sequences analyzed ..: " . sprintf("%${strLen}s\n", commify($readCount+$readsUnmp));
$summary .= "    Mapped Reads ................: " . sprintf("%${strLen}s (%.1f%%)\n", commify($readCount), $frac{validReads});
$summary .= "    Marked \"Duplicated\" .........: " . sprintf("%${strLen}s (%.1f%%)\n", commify($dupReads), $frac{dupReads});
$summary .= "    Mapped Read Pairs ...........: " . sprintf("%${strLen}s (%.1f%%)\n", commify($mappedPairs), $frac{mappedPairs}*2);
$summary .= "    Unmapped Reads ..............: " . sprintf("%${strLen}s (%.1f%%)\n", commify($readsUnmp), $frac{unmappedReads});
#$summary .= "  Mapped Reads -> Bases .........: " . sprintf("%${strLen}s (%s)\n", commify($totalLength), sizeBases($totalLength));
$summary .= "  Average Read Length ...........: " . sprintf("%${strLen}.0f\n", $totalLength/$readCount);
$summary .= "\n\n";
$summary .= "Alignment Statistics\n";
$summary .= "--------------------\n";
$summary .= "\n";
$summary .= "  (M) Aligned Bases .............: " . sprintf("%${strLen}s (%.1f%%)\n", commify($cigar{M}), $frac{M});
$summary .= "  (S) Soft Clipped Bases ........: " . sprintf("%${strLen}s (%.1f%%)\n", commify($cigar{S}), $frac{S});
$summary .= "  (I) Insertions ................: " . sprintf("%${strLen}s (%.1f%%)\n", commify($cigar{I}), $frac{I});
$summary .= "  (D) Deletions .................: " . sprintf("%${strLen}s (%.1f%%)\n", commify($cigar{D}), $frac{D});
$summary .= "  (N) Skipped Bases .............: " . sprintf("%${strLen}s (%.1f%%)\n", commify($cigar{N}), $frac{N});
$summary .= "  (H) Hard Clipped Bases ........: " . sprintf("%${strLen}s (%.1f%%)\n", commify($cigar{H}), $frac{H});
$summary .= "  (P) Padded Bases ..............: " . sprintf("%${strLen}s (%.1f%%)\n", commify($cigar{P}), $frac{P});
$summary .= "  (=) Sequence Matches ..........: " . sprintf("%${strLen}s (%.1f%%)\n", commify($cigar{'='}), $frac{'='});
$summary .= "  (X) Sequence Mismatch .........: " . sprintf("%${strLen}s (%.1f%%)\n", commify($cigar{X}), $frac{X});
$summary .= "\n\n";
$summary .= "Clipping Statistics\n";
$summary .= "-------------------\n";
$summary .= "\n";

if($clip5p) {

   $modalLength = (defined $stat5p->mode()) ? $stat5p->mode() : "No count differences.";

   $summary .= "  Clipped Bases at 5' End .......: " . sprintf("%${strLen}s (%.1f%%)\n", commify($clip5p), $frac{clipped5p});
   $summary .= "   Clipped Length Stats >-\n";
   $summary .= "      Average ............: " . sprintf("%4.0f\n",$stat5p->mean());
   $summary .= "      Median .............: " . sprintf("%4.0f\n",$stat5p->median());
   $summary .= "      StdDev .............: " . sprintf("%4.0f\n",$stat5p->standard_deviation());
   $summary .= "      Modal ..............: " . sprintf("%4.0f\n",$modalLength);
   $summary .= "      Min ................: " . sprintf("%4.0f\n",$stat5p->min());
   $summary .= "      Max ................: " . sprintf("%4.0f\n",$stat5p->max());
   $summary .= "\n";
}
else {
   $summary .= "No 5' soft clipping applied.\n";
}

if($clip3p) {

   $modalLength = (defined $stat3p->mode()) ? $stat3p->mode() : "No count differences.";

   $summary .= "  Clipped Bases at 3' End .......: " . sprintf("%${strLen}s (%.1f%%)\n", commify($clip3p), $frac{clipped3p});
   $summary .= "   Clipped Length Stats ->\n";
   $summary .= "      Average ............: " . sprintf("%4.0f\n",$stat3p->mean());
   $summary .= "      Median .............: " . sprintf("%4.0f\n",$stat3p->median());
   $summary .= "      StdDev .............: " . sprintf("%4.0f\n",$stat3p->standard_deviation());
   $summary .= "      Modal ..............: " . sprintf("%4.0f\n",$modalLength);
   $summary .= "      Min ................: " . sprintf("%4.0f\n",$stat3p->min());
   $summary .= "      Max ................: " . sprintf("%4.0f\n",$stat3p->max());
   $summary .= "\n";
}
else {
   $summary .= "No 3' soft clipping applied.\n";
}

## finally finished
$timeEnd = time();

$summary .= "\n";
$summary .= "Processing Statistics\n";
$summary .= "---------------------\n";
$summary .= "\n";
$summary .= "  Alignment File Type ....: " . $config{alnType} . "\n";
$summary .= "  Processing Rate ........: " . commify($processingRate) . " reads per second\n";
$summary .= "  Total Time .............: " . duration($timeEnd-$timeStart) . "\n";
$summary .= "\n";


print $summary;



### ---------------------------------------------------
sub commify {
### ---------------------------------------------------
   my $number = shift;
   $number = reverse $number;
   $number =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
   $number = reverse $number;
   return $number;
}

### ---------------------------------------------------
sub readsPerSecond {
### ---------------------------------------------------
   my ($noOfReads,$timeToProcess)=@_;
   $timeToProcess = ($timeToProcess<1) ? return 'ALL' : return (int($noOfReads/$timeToProcess));
}

### ---------------------------------------------------
sub sizeBases {
### ---------------------------------------------------
   my $bases=shift;
   return sprintf("%.1f Gb",$bases/1_000_000_000);
}

### ---------------------------------------------------
sub duration {
### ---------------------------------------------------
   my $total=shift;

   ## it usually takes minutes
   my $min=int($total/60);
   my $sec=$total-($min*60);

   return $min==0 ? "$total sec" : "$min min $sec sec";
}

### ---------------------------------------------------
sub calculateProcessRate {
### ---------------------------------------------------
   my ($readCount,$timeStart) = @_;
   my $timeEnd = time();

   if( ($timeEnd-$timeStart)<1 ) {
      return 0;
   }
   else {
      return (int($readCount/($timeEnd-$timeStart)));
   }
}


__DATA__

SAM File :
-rw-r----- 1 klages abt_srv 171G 2014.10.15 14:19:25 core_L4637-1_2Leber-A_NxSeq_marMar2.sam

Input Statistics
----------------

  Number of sequences analyzed ..: 446,743,644
       Valid CIGAR  .............: 410,453,052 (91.9%)
         Total Number of Bases ..: 61,157,504,748
       Skipped Sequences (*) ....: 36,290,592 (8.1%)


Alignment Statistics
--------------------

  (M) Aligned Bases .............: 46,352,306,693 (75.8%)
  (S) Soft Clipped Bases ........: 14,800,026,203 (24.2%)
  (I) Insertions ................: 5,171,852 (0.0%)
  (D) Deletions .................: 6,828,217 (0.0%)
  (N) Skipped Bases .............: 0 (0.0%)
  (H) Hard Clipped Bases ........: 0 (0.0%)
  (P) Padded Bases ..............: 0 (0.0%)
  (=) Sequence Matches ..........: 0 (0.0%)
  (X) Sequence Mismatch .........: 0 (0.0%)


Clipping Statistics
-------------------

  Clipped Bases at 5' End : 7,329,469,546 (12.0%)
   Clipped Length Stats :
      Average ....: 58
      Median .....: 50
      StdDev .....: 48
      Modal ......: 1
      Min ........: 1
      Max ........: 148

  Clipped Bases at 3' End :  7,470,556,657 (12.2%)
   Clipped Length Stats:
      Average ....: 58
      Median .....: 50
      StdDev .....: 48
      Modal ......: 1
      Min ........: 1
      Max ........: 148


real  158m16.714s
user  150m40.436s
sys   18m41.717s
