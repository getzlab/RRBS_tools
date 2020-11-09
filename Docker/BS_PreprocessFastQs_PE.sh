if [ $# -eq 0 ]; then
    echo "Usage: $0 -1 FQ1 -2 FQ2 -s SAMPLE -g GENOME -r REFERENCE_FA -z CHROM_SIZES [-c C_THREADS -n -a -o]"
    echo ""
    echo "FQ1       : first fastq file of a pair"
    echo "FQ2       : second fastq file of a pair"
    echo "SAMPLE    : sample name, used for output naming"
    echo "GENOME    : e.g. 'mm9'/'mm10' (mouse) or 'hg19'/'hg38' (human)."
    echo "C_THREADS : computing(alignment) threads to use."
    echo "m         : total memory (gigs) "
    echo "n         : special flag for 2-color devices NextSeq and NovaSeq."
    echo "o         : overlap-based clipping using bamutils (experimental)"
    #echo "r         : Reference FASTA (Must not be compressed)"
    echo ""
    echo " Version : $__VERSION"
    exit
fi

## parameters
while getopts 1:2:i:s:g:c:r:q:z:m:nx option
do
    case "${option}"
    in
	1) FQ1_INPUT=${OPTARG};;
	2) FQ2_INPUT=${OPTARG};;
	s) SAMPLE=${OPTARG};;
	g) GENOME=${OPTARG};;
	c) C_THREADS=${OPTARG};;
	m) mem_gigs=${OPTARG};;
	n) IS_TWO_COLOR_SEQ=1;;
	#r) REFERENCE_FA=${OPTARG};;
	x) CUTADAPT=1;;
    esac
done

#########################
##STEP 1
##inputs
#-1 $FQ1
#-2 $FQ2
#-g genome
#-r genome_reference_fa
#-s $samplename
#-c $threads (task-specific, 4 probably ideal)
#-m $mem_gigs (task-specific)
#-x cutadapt (optional)
#-n istwo color seq (optional)


function print_str {
    str=$1
    echo ""
    echo $(date +"### [%Y-%m-%d %H:%M:%S] $str")
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
#if [ -z ${genome_index[$GENOME]+"is_set_value"} ];then
#    echo "[ERROR] Unknown reference genome '$GENOME'."
#    exit
#fifi

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

if [ $C_THREADS -lt 4 ]; then
    echo "[ERROR] Use 4 or more threads for computing, please."
    exit
fi

## PARAMETER CHECK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## avoid funny locales for core components QC, mapping, transcript assembly
export LC_ALL="POSIX"




## START WORKFLOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ----------------------------------------------------------------
## (1) Prepare working/results directory
## ----------------------------------------------------------------
results_dir=${working_dir}/${SAMPLE}_${GENOME}_PE
fastqc_results_dir=$results_dir/00_FastQC
mkdir -p $fastqc_results_dir
cd $results_dir


## ----------------------------------------------------------------
## (2) pre-QC
## ----------------------------------------------------------------


## PARAMETER CHECK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## (1,2) Fastq input data
if [ ! -r "$FQ1_INPUT" ] || [ ! -r "$FQ2_INPUT" ]; then
    echo "$FQ1_INPUT or $FQ2_INPUT cannot be read."
    exit
fi

## absolute path
FQ1=$(realpath -s $FQ1_INPUT)
FQ2=$(realpath -s $FQ2_INPUT)


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

if [ $C_THREADS -lt 4 ]; then
    echo "[ERROR] Use 4 or more threads for computing, please."
    exit
fi

## PARAMETER CHECK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## avoid funny locales for core components QC, mapping, transcript assembly
export LC_ALL="POSIX"


## START WORKFLOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ----------------------------------------------------------------
## (1) Prepare working/results directory
## ----------------------------------------------------------------
results_dir=${working_dir}/${SAMPLE}_${GENOME}_PE
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

if [ $CUTADAPT  ]; then
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
    echo User specified to cut adaptor sequences
    echo current adaptor sequences: ${adaptors[illumina]}
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

fi
