#!/bin/bash

#adunford:making some essential edits to get this to work in docker/google cloud environment
#namely changing the way it references the reference genome information, which currently is hardcoded as if it is in a static local area
#and changing the binaries to point to the default path locations installed in the docker instead of the hardcoded areas listed

:<<'WGBS_WORKFLOW'

  * PE Data only

  This workflow works with libraries made with Illumina swift kit which adds a random sequence "low complexity tail"
  during adapter ligation see here: https://swiftbiosci.com/wp-content/uploads/2016/09/16-0853-Tail-Trim-TN.pdf
  For these it is necessary to trim 10bp from the start of each read.

  ALIGNER: bsmap
  CALLER : mcall

WGBS_WORKFLOW

__VERSION='1.0 - 11/2018'

## default values
C_THREADS=40 # compute threads
p_threads=4  # respect samtools pipe

## mainly for GATK
ulimit -n 4096

if [ $# -eq 0 ]; then
  echo "Usage: $0 -1 FQ1 -2 FQ2 -s SAMPLE -g GENOME -r REFERENCE_FA -z SIZE [-c C_THREADS -n -a -o]"
  echo ""
  echo "FQ1       : first fastq file of a pair"
  echo "FQ2       : second fastq file of a pair"
  echo "SAMPLE    : sample name, used for output naming"
  echo "GENOME    : either 'mm9'/'mm10' (mouse) or 'hg19'/'hg38' (human)."
  echo "C_THREADS : computing(alignment) threads to use."
  echo "n         : special flag for 2-color devices NextSeq and NovaSeq."
  echo "a         : perform not only CpG, but additionally CpA,C,T calling"
  echo "o         : overlap-based clipping using bamutils (experimental)"
  echo ""
  echo " Version : $__VERSION"
  exit
fi

## parameters
while getopts 1:2:i:s:g:c:r:z:nao option
do
    case "${option}"
        in
        1) FQ1_INPUT=${OPTARG};;
        2) FQ2_INPUT=${OPTARG};;
        s) SAMPLE=${OPTARG};;
        g) GENOME=${OPTARG};;
        c) C_THREADS=${OPTARG};;
        n) IS_TWO_COLOR_SEQ=1;;
	r) REFERENCE_FA=${OPTARG};;
	z) SIZE=${OPTARG};;
        a) DO_ALL=1;;
        o) DO_OBC=1;;
    esac
done

## safer way
set -e

## common stuff
working_dir=$(pwd)
date_flag=$(date +%Y%m%d)  # 20181025.144207
temp_dir=temp.WGBS.$date_flag
declare -a bio_muelltonne  # array for final cleanup

## Effective Genome Sizes (for bamCoverage)
## https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
declare -A egs
egs["hg19"]=2864785220
egs["hg38"]=2913022398
egs["mm9"]=2620345972
egs["mm10"]=2652783500


## Genome Name Lookup
#adunford: Instead of making hard-coded environment variables just require user too provide FASTA
#declare -A genome_index
#genome_index[hg19]='/project/genomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa'
#genome_index[hg38]='/project/genomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'
# genome_index[mm9]='/project/genomes/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa'
#genome_index[mm10]='/project/genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa'
genome_index=${REFERENCE_FA}


## Genome Size Lookup
#adunford: as above.
#declare -A chrom_sizes
#chrom_sizes[hg19]='/project/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.chrom.sizes'
#chrom_sizes[hg38]='/project/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.chrom.sizes'
# chrom_sizes[mm9]='/project/genomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.chrom.sizes'
#chrom_sizes[mm10]='/project/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.chrom.sizes'
chrom_sizes=${SIZE}

#adunford: might need to expose this as an input parameter
## Adaptors to trim ()
declare -A adaptors
adaptors[illumina]='AGATCGGAAGAGC'
adaptors[smallrna]='TGGAATTCTCGG'
adaptors[nextera]='CTGTCTCTTATA'

## be verbose
function print_str {
  str=$1
  echo ""
  echo $(date +"### [%Y-%m-%d %H:%M:%S] $str")
}

function trim_galore {
    echo "TrimGalore clipping incl RRBS alternative to Nugen's scripts - may run multi-threaded."
}

function bed2bigbed {
    bed_in=$1
    bed_bb=${bed_in%.bed}.bb
    b2bb_bin=/package/sequencer/bin/bedToBigBed
    reformat=/package/sequencer/bin/mcall_bed_reformat.pl

    ## This is incredibly slow :-(
    #    ## temporary file
    #    tmp_file='.tmp_b2bb'
    #    rm -f $tmp_file
    #
    #    ## mainly for number format
    #    export LC_ALL=POSIX
    #
    #    # IFS=$'\t' : word separator, special with $ is syntax for string literals with escape sequences
    #    # -r        : do not interprete escape characters like '\'
    #    # -a        : read line and split words into array 'line'
    #    while IFS=$'\t' read -r -a line; do
    #        # skip header
    #        [[ ${line[0]} =~ ^# ]] && continue
    #
    #        # totalC > 4
    #        [[ ${line[4]} -lt 5 ]] && continue
    #
    #        # floating point arithmetic not supported by bash, using 'bc'
    #        # convert result to integer
    #        printf -v frac "%.0f" $(echo "${line[3]}*100" | bc)
    #
    #        # output columns (format) have been taken from original "EPP",
    #        # pipeline / Michael Ziller. mcall bed output is sorted.
    #        printf "%s\t%u\t%u\t%s\t%u\t\t+\n" \
    #          ${line[0]} \
    #          ${line[1]} \
    #          ${line[2]} \
    #          "'${frac}%[${line[3]}]'" \
    #          $frac \
    #          >> $tmp_file
    #
    #    done < $bed_in

    tmp_file=$($reformat $bed_in)

    ## convert BED to BigBED using bedToBigBed
    $b2bb_bin \
      $tmp_file \
      ${chrom_sizes[$GENOME]} \
      $bed_bb

    ## remove tmp file here
    rm -f $tmp_file
}


## PARAMETER CHECK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## (1,2) Fastq input data
if [ ! -r "$FQ1_INPUT" ] || [ ! -r "$FQ2_INPUT" ]; then
    echo "$FQ1_INPUT or $FQ2_INPUT cannot be read."
    exit
fi

## absolute path
FQ1=$(realpath -s $FQ1_INPUT)
FQ2=$(realpath -s $FQ2_INPUT)



## (g) Genome Reference
## returns empty string if $genome_index[$GENOME] is NOT set, "is_set_value" if set.
## --> bash parameter substitution
if [ -z ${genome_index[$GENOME]+"is_set_value"} ];then
    echo "[ERROR] Unknown reference genome '$GENOME'."
    exit
fi

if [ ! -r ${chrom_sizes[$GENOME]} ];then
    echo "[ERROR] '${chrom_sizes[$GENOME]}' not found for reference genome '$GENOME'."
    exit
fi

## (s) Sample Name
if [ -z $SAMPLE ]; then
    echo "[ERROR] Please provide SAMPLE name."
    exit
fi

## (c) Threads - some basic testing for positive integer
if [[ ! $C_THREADS =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] Not a valid number '$C_THREADS' for computing threads."
    exit
fi

if [ $C_THREADS -lt 5 ]; then
    echo "[ERROR] Use 5 or more threads for computing, please."
    exit
fi

## PARAMETER CHECK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## avoid funny locales for core components QC, mapping, transcript assembly
export LC_ALL="POSIX"


## programs to use
#       pigz=/package/sequencer/bin/pigz
#      gatk4=/package/sequencer/bin/gatk
#      bsmap=/package/sequencer/bin/bsmap
#      mcall=/package/sequencer/bin/mcall
#     fastqc=/package/sequencer/bin/fastqc
#    multiqc=/package/sequencer/bin/multiqc
#   bamstats=/package/sequencer/bin/bamStat.pl
#   bamutils=/package/sequencer/bin/bamutils
#   cutadapt=/package/sequencer/bin/cutadapt
#   samtools=/package/sequencer/bin/samtools
#bamCoverage=/package/sequencer/bin/bamCoverage

## START WORKFLOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ----------------------------------------------------------------
## (1) Prepare working/results directory
## ----------------------------------------------------------------
results_dir=${working_dir}/wgbs_${SAMPLE}_${GENOME}_PE_${date_flag}
fastqc_results_dir=$results_dir/00_FastQC
mkdir -p $fastqc_results_dir
cd $results_dir


## ----------------------------------------------------------------
## (2) pre-QC
## ----------------------------------------------------------------
print_str "FastQC : pre-QC"

fastqc \
  --outdir $fastqc_results_dir \
  --threads $C_THREADS \
  --dir $fastqc_results_dir \
  $FQ1 $FQ2


echo "[INFO] Applying Adaptor and Quality Trimming to Input Data"

## ----------------------------------------------------------------
## (3) Adaptor and Quality Trimming
## ----------------------------------------------------------------
print_str "Running cutadapt - Trimming Illumina adaptors and Quality"

min_length=25
min_qual=20
swift_trim_off=10

# output names, no absolute paths necessary
fq1_trimmed=${SAMPLE}_R1.trm.fq.gz
fq2_trimmed=${SAMPLE}_R2.trm.fq.gz

# If data has been generated on NextSeq or NovaSeq 2-color machines,
# we need to apply a different type of quality clipping (due to dark cycles)
trimming_type='--quality-cutoff'
if [ $IS_TWO_COLOR_SEQ ] ; then
    echo "[INFO] Applying 2-color-device trimming mode"
    trimming_type='--nextseq-trim'
fi

## R1----->
## ~~~~~~~|=======================|~~~~~~~
##   AD   |        DNA            |   AD
## ~~~~~~~|=======================|~~~~~~~
##                                <-----R2

# quality clipping is done BEFORE adaptor clipping
cutadapt \
  $trimming_type $min_qual \
  --overlap 5 \
  --minimum-length $(($min_length+2*$swift_trim_off)) \
  --cores $(($C_THREADS/2)) \
  --adapter ${adaptors[illumina]} \
  -A ${adaptors[illumina]} \
  --interleaved \
  $FQ1 \
  $FQ2 \
  2>.log_a | \
cutadapt \
  --interleaved \
  --minimum-length $min_length \
  --cores $(($C_THREADS/2)) \
  --cut $swift_trim_off \
  --cut -$swift_trim_off \
  -U $swift_trim_off -U -$swift_trim_off \
  --output $fq1_trimmed \
  --paired-output $fq2_trimmed \
  - \
  &>.log_b

## merge log output
cat .log_a .log_b > ${SAMPLE}_FastQ-Trimming-Adp.log

bio_muelltonne+=($fq1_trimmed $fq2_trimmed .log_a .log_b)


## ----------------------------------------------------------------
## (4) post-QC
## ----------------------------------------------------------------
print_str "FastQC : post-QC"

fastqc \
  --outdir $fastqc_results_dir \
  --threads $C_THREADS \
  --dir $fastqc_results_dir \
  $fq1_trimmed $fq2_trimmed

## finally set FQ1
FQ1=$fq1_trimmed
FQ2=$fq2_trimmed

pwd
echo $FQ1 $FQ2
## ----------------------------------------------------------------
## (5) Alignment - bsmap
## ----------------------------------------------------------------
print_str "bsmap Alignment : $SAMPLE on $GENOME (${genome_index[$GENOME]})"

bam_file=${SAMPLE}_${GENOME}.bsmap.srt.bam

# -s = seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.
# In pipeline this was set to 1 (WGBS mode); here I revert it back to default, 12
#
# -u = report unmapped reads, default=off
# -R = print corresponding reference sequences in SAM output, default=off
#
seed_size=16 # default=16(WGBS mode), 12(RRBS mode). min=8, max=16
max_ins=1000 # max insert size for PE mapping (-x)

echo bsmap  -v 0.1 -s $seed_size -q 20 -w 100 -S 1 -u -R -x $max_ins -p $(($C_THREADS-$p_threads)) -d ${genome_index} -a $FQ1  -b $FQ2 samtools sort -m 8G --threads $p_threads --output-fmt BAM -o $bam_file -

#need to make this modular
mem_gigs=2

bsmap \
  -v 0.1 -s $seed_size -q 20 -w 100 -S 1 -u -R \
  -x $max_ins \
  -p $(($C_THREADS-$p_threads)) \
  -d ${genome_index} \
  -a $FQ1 \
  -b $FQ2 \
  | samtools sort -m ${mem_gigs}G --threads $p_threads --output-fmt BAM -o $bam_file -

print_str "Creating BAM Index for $bam_file"
samtools index -@ $((2*$p_threads)) $bam_file

print_str "Creating BAM Statistcis $bam_file"
/src/bamStat.pl --bam $bam_file --threads $((2*$p_threads)) --verbose | tee ${bam_file}.stats.txt

bio_muelltonne+=($bam_file ${bam_file}.bai)


## ----------------------------------------------------------------
## (6) Overlap-based Clipping -- experimental !
## ----------------------------------------------------------------
if [ $DO_OBC ]; then

    bam_file_oc=${bam_file%.bam}.oc.bam

    bamutils \
      clipOverlap \
      --in $bam_file \
      --out $bam_file_oc \
      --stats \
      --params \
      --unmapped \
      --noPhoneHome \
      --poolSize 10000000 # def 1000000

    print_str "Creating BAM Index for $bam_file_oc"
    samtools index -@ $((2*$p_threads)) $bam_file_oc

    print_str "Creating BAM Statistcis $bam_file_oc"
    bamstats --bam $bam_file_oc --threads $((2*$p_threads)) --verbose | tee ${bam_file_oc}.stats.txt

    bio_muelltonne+=($bam_file_oc ${bam_file_oc}.bai)

    # set bam_file to new one
    bam_file=$bam_file_oc
fi


## ----------------------------------------------------------------
## (7) Deduplication
## ----------------------------------------------------------------
max_open_files=$(ulimit -n) # $max_open_files/1000*1000 = 4096/1000*1000 = 4000 :-) ++
mhs_mem=16g # max heap size
bam_dedup=${bam_file%.bam}.rd.bam

## We skip optical duplicate finding by setting READ_NAME_REGEX=null.
## With HiSeq 4000 and NovaSeq 6000 there shouldn't be any "optical duplicates"
## due to their patterned flowcells. We will run only few WGBS libs on other machines.
## sklages, 2018-11-15

gatk \
  MarkDuplicates \
  --java-options "-XX:ParallelGCThreads=$C_THREADS -Xmx$mhs_mem -Djava.io.tmpdir=$temp_dir" \
  --ASSUME_SORT_ORDER=coordinate \
  --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$(($max_open_files/1000*1000)) \
  --REMOVE_DUPLICATES=true \
  --MAX_RECORDS_IN_RAM=10000000 \
  --INPUT=$bam_file \
  --OUTPUT=$bam_dedup \
  --METRICS_FILE=${bam_dedup}.dedup-metrics.txt \
  --TMP_DIR=$temp_dir \
  --VALIDATION_STRINGENCY=LENIENT \
  --READ_NAME_REGEX=null

print_str "Creating BAM Index for $bam_dedup"
samtools index -@ $((2*$p_threads)) $bam_dedup

print_str "Creating BAM Statistcis $bam_dedup"
/src/bamStat.pl --bam $bam_dedup --threads $((2*$p_threads)) --verbose | tee ${bam_dedup}.stats.txt

# set bam_file to new one - keep final BAM
bam_file=$bam_dedup


## ----------------------------------------------------------------
## (8) Methylation Calling
## TODO: maybe compress BED output; files may get huge for CpATC
## ----------------------------------------------------------------
print_str "Running Methylation Calling CpG"

prefix=${SAMPLE}_${GENOME}
mcall \
  --threads $C_THREADS \
  --reference ${genome_index[$GENOME]} \
  --sampleName $prefix \
  --mappedFiles $bam_file \
  --outputDir $results_dir \
  --webOutputDir $results_dir

## output generated by mcall
## run.config.64964
## AM-RRBS-88_mm9.G.bed -> AM-WGBS-88_mm9.bsmap.srt.oc.rd.bam.G.bed
## AM-WGBS-88_mm9.bsmap.srt.oc.rd.bam.G.bed
## AM-WGBS-88_mm9.bsmap.srt.oc.rd.bam.HG.bed
## AM-WGBS-88_mm9.bsmap.srt.oc.rd.bam_stat.txt
## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
## <---        $bam_file         --->
cpx_bed=${bam_file}.CpG.bed
mv ${bam_file}.G.bed $cpx_bed

cpx_txt=${bam_file}.CpX.stats.txt
mv ${bam_file}_stat.txt $cpx_txt

## converting bed to bb
print_str "Converting BED ($cpx_bed) to BigBed"
bed2bigbed $cpx_bed

bio_muelltonne+=(${prefix}.G.bed ${bam_file}.HG.bed)

## in case someone wants to have non-CpG calls
if [ $DO_ALL ] ; then
    print_str "Running Methylation Calling on non-CpG bases"

    for base in A C T;do
        rm -f ${prefix}.G.bed # ln: failed to create symbolic link 'AM-RRBS-88_mm9.G.bed': File exists
        mcall \
          --threads $C_THREADS \
          --reference ${genome_index} \
          --sampleName $prefix \
          --mappedFiles $bam_file \
          --outputDir $results_dir \
          --webOutputDir $results_dir \
          --reportCpX $base

        ## rename file
        cpx_bed=${bam_file}.Cp${base}.bed
        mv ${bam_file}.G.bed $cpx_bed

        ## this is some common stats, all the same for all CpACGT
        ## so we simply overwrite :-)
        cpx_txt=${bam_file}.CpX.stats.txt
        mv ${bam_file}_stat.txt $cpx_txt

        ## converting bed to bb
        print_str "Converting BED ($cpx_bed) to BigBed"
        bed2bigbed $cpx_bed
    done
fi


## [..] abort further execution because Python 3 was configured to use ASCII as encoding for the environment.
## This system lists a couple of UTF-8 supporting locales that you can pick from.
## deeptools,multiqc
export LC_ALL=de_DE.UTF-8

## ----------------------------------------------------------------
## (10) Create MultiQC Report
## ----------------------------------------------------------------
print_str "Running MultiQC"

multiqc \
  --title ${SAMPLE}_${GENOME} \
  $results_dir


## ----------------------------------------------------------------
## (11) Cleaning up
## ----------------------------------------------------------------
print_str "Cleaning up"

rm -f run.config.*
bio_muelltonne+=($temp_dir)

for item in "${bio_muelltonne[@]}"
do
    echo " removing '$item'"
    rm -fr $item
done

echo ""
echo "***********************************************************"
echo "***     Problem:"
echo "***     https://github.com/marcelm/cutadapt/issues/366"
echo "***********************************************************"
echo ""
