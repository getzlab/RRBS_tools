###??? To do
# run fastqc pre-trimmed, post-trimmed (parallel to bsmap)
# preprocess_fastqs
# -- May want threads-1 (this seems to have been the setting)
# Cutadapt:
# --- Understand swift_trim_off=10 (are we trimming too much?)
# --- Check if our data was on "NextSeq or NovaSeq 2-color machines" (if so needs different quality trimming)

# mcall
# -- understand bed vs bigbed
# -- decide if to gzip beds (likely yes!)

# docker not optional in tasks, only in workflow; hardcode version

task fastqc{
    File fastq1
    File fastq2
    String? fastq_suffix='.fq.gz'
    String? prefix1=basename(fastq1, fastq_suffix)
    String? prefix2=basename(fastq2, fastq_suffix)

    String? outdir_suffix=""
    String? fastqc_args=""

    Int threads
    Int seed_size #seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.
    String? docker="adunford/bsmap_to_mcall:latest"
    Int disk_size_gb
    Int mem
    Int preemtible

    command {
        outdir="fastqc"
        if [ -n "${outdir_suffix}" ]; then
            outdir="$outdir_${outdir_suffix}"
        fi
        mkdir -p $outdir
        $fastqc \
          --outdir $outdir \
          --threads ${threads} \
          ${fastqc_args} \
          ${fastq1} ${fastq2}
    }

    runtime {
    	    docker: "${docker}"
    	    memory: "${mem}GB"
            cpu: "${threads}+1"
    	    disks: "local-disk ${disk_size_gb} HDD"
    	    preemptible: "${preemtible}"
    }
    output {
        # fastqc outdir.tar.gz ??? (to send to multiqc)
        #CRC-0021-T-00_5_R1.trm_fastqc.html
        #CRC-0021-T-00_5_R2.trm_fastqc.zip
    }
}

task preprocess_fastqs{
    File fastq1
    File fastq2
    String sample_id
    #String genome_reference #@adunford, this is not needed in Terra (all defined in config)
    File reference_fa
    File reference_sizes #chromosome sizes

    Int threads
    Int seed_size #seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.
    String docker
    Int disk_size_gb
    Int mem
    Int preemtible

    command {
    	    /src/NEW_WORKFLOW_NAME-???.sh -1 ${fastq1} -2 ${fastq2} -s ${sample_id} -c ${threads}

            ### cutadapt
            ### cutadapt
 	    }

    runtime {
    	    docker: "${docker}"
    	    memory: "${mem}GB"
            cpu: "${threads}"
    	    disks: "local-disk ${disk_size_gb} HDD"
    	    preemptible: "${preemtible}"
    }
    output {
	    File fastq1_trimmed = "${sample_id}_R1.trm.fq.gz"
	    File fastq2_trimmed = "${sample_id}_R2.trm.fq.gz"
    }
}

task bsmap{
    File fastq1
    File fastq2
    String sample_id
    #GENOME    : either 'mm9'/'mm10' (mouse) or 'hg19'/'hg38' (human).
    String genome_reference
    File reference_fa
    File reference_sizes
    Int? seed_size="12" #default=12(RRBS mode), 16(WGBS mode). min=8, max=16.
    Int? max_insert_size="1000" # max insert size for PE mapping (-x)
    String? bsmap_args="-q 20 -w 100 -S 1 -u -R" #??? understand each param in context of RRBS vs WGBS

    Int bsmap_mem
    Int bsmap_threads="12"
    Int? samtools_pipe_threads="4"

    String docker
    Int mem
    Int threads
    Int disk_size_gb
    Int preemtible

    command {
            # -s = seed size, default=12(RRBS mode), 16(WGBS mode). min=8, max=16
            # -u = report unmapped reads, default=off
            # -R = print corresponding reference sequences in SAM output, default=off
            # -q = #???
            # -w = #???
            # -S = 1 #@adunford ??? - stranded? Verify good for all data types
            echo total mem: ${mem}
            echo mem per samtools thread: $(((${mem}-1)/${samtools_pipe_threads}))G

            bam_file=${sample_id}.bsmap.srt.bam

            ###@adunford, please review the memory and thread allocations

            bsmap \
              -v 0.1 -s $seed_size ${bsmap_args} \
              -x ${max_insert_size} \
              -p $((${threads}-$samtools_pipe_threads)) \
              -d ${genome_index} \
              -a $FQ1 \
              -b $FQ2 \
              | samtools sort -m  $(((${mem}-${bsmap_mem}-1)/$samtools_pipe_threads))G --threads $samtools_pipe_threads --output-fmt BAM -o $bam_file

            print_str "Creating BAM Index for $bam_file"
            samtools index -@ ${threads} $bam_file

            print_str "Creating BAM Statistcis $bam_file"
            /src/bamStat.pl --bam $bam_file --threads $((${threads}-1)) --verbose | tee ${sample_id}.bsmap.srt.stats.txt
 	    }

    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${preemtible}"
    }
    output {
        File bam = "${sample_id}.bsmap.srt.bam"
        File bam_index = "${sample_id}.bsmap.srt.bam.bai"
        File bam_stats = "${sample_id}.bsmap.srt.stats.txt"
    }
}

task clip_overlaps {
    File bam_file
    File bam_index
    String? prefix=basename(bam_file, ".bam")

    String? clip_overlap_args="--stats --params --unmapped --noPhoneHome --poolSize 10000000"

    String docker
    Int mem
    Int threads
    Int disk_buffer_gb = "20"
    Int disk_size_gb = ceil(size(bam_file, "G") + size(bam_index, "G")) * 3 #for input, intermediates, output
    Int preemtible

    command {
        bam_file_oc=${bam_file%.bam}.oc.bam

        bamutils clipOverlap --in ${bam_file} --out $bam_file_oc ${clip_overlap_args}

        echo $(date +"### [%Y-%m-%d %H:%M:%S] Creating BAM Index for $bam_file_oc")
        samtools index -@ ${threads} $bam_file_oc

        echo $(date +"### [%Y-%m-%d %H:%M:%S] Creating BAM Statistics $bam_file_oc")
        bamstats --bam $bam_file_oc --threads $(${threads}-1) --verbose | tee ${bam_file_oc}.stats.txt
    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${preemtible}"
    }

    output {
        bam_overlap_clipped="${prefix}.oc.bam" #@adunford: does this work?
        bam_overlap_clipped_index="${prefix}.oc.bam.bai"
        bam_overlap_clipped_stats="${prefix}.oc.bam.stats.txt"
    }
}

task markduplicates {
    File bam_file
    File bam_index
    String sample_id
    String? prefix=basename(bam_file, ".bam")

    String remove_dups="true" #WARNING: "true" LEADS TO LOSS OF READS IN FINAL BAM
    Int? max_records_in_ram="500000" #As in GTEx; Helene: "10000000"
    Int? max_open_files="8000" #GATK default
    Int? java_mem="3" # max heap size #As in GTEx; Helene: "16G"
    String? read_name_regex="null" #skips optical duplicate finding (this is setting for HiSeq 4000 and NovaSeq 6000 that shouldnt have optical duplicates)

    String docker
    Int mem
    Int threads
    Int disk_buffer_gb = "20"
    Int disk_size_gb = ceil(size(bam_file, "G") + size(bam_index, "G")) * 3 #for input, intermediates, output
    Int preemtible

    command {
        bam_dedup=${bam_file%.bam}.md.bam

        ## We skip optical duplicate finding by setting READ_NAME_REGEX=null.
        ## With HiSeq 4000 and NovaSeq 6000 there shouldn't be any "optical duplicates"
        ## due to their patterned flowcells. We will run only few WGBS libs on other machines.
        ## sklages, 2018-11-15
        temp_dir="tempdir"
        mkdir $temp_dir

        gatk \
          MarkDuplicates \
          --java-options "-XX:ParallelGCThreads=${threads} -Xmx${java_mem}g -Djava.io.tmpdir=$temp_dir" \
          --ASSUME_SORT_ORDER=coordinate \
          --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=${max_open_files} \
          --REMOVE_DUPLICATES=${remove_dups} \
          --MAX_RECORDS_IN_RAM=${max_records_in_ram} \
          --INPUT=${bam_file} \
          --OUTPUT=$bam_dedup \
          --METRICS_FILE=${bam_dedup}.dedup-metrics.txt \
          --TMP_DIR=$temp_dir \
          --VALIDATION_STRINGENCY=LENIENT \
          --READ_NAME_REGEX=${read_name_regex}

        echo $(date +"### [%Y-%m-%d %H:%M:%S] Creating BAM Index for $bam_dedup")
        samtools index -@ ${threads} $bam_dedup

        echo $(date +"### [%Y-%m-%d %H:%M:%S] Creating BAM Statistics $bam_dedup")
        /src/bamStat.pl --bam $bam_dedup --threads $(${threads}-1) --verbose | tee ${bam_dedup}.stats.txt #BK: allow 1 thread for "tee"
    }

    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${preemtible}"
    }
    output {
        bam_md="${prefix}.md.bam"
        bam_md_index="${prefix}.md.bam.bai"
        bam_md_stats="${prefix}.md.bam.stats.txt"
    }
}

task mcall {
    File bam_file
INPUTS ???

    command {
        echo $(date +"### [%Y-%m-%d %H:%M:%S] Running Methylation Calling CpG")

        prefix=${sample_id}

        mcall \
          --threads ${threads} \
          --reference ${genome_index[$GENOME]} \
          --sampleName ${sample_id} \ #BK: removed genome index
          --mappedFiles ${bam_file} \
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
        cpg_bed=${bam_file}.CpG.bed
        mv ${bam_file}.G.bed $cpg_bed

        #Stats are not unique to CpG (would be the same for e.g. CpA analysis), hence naming CpX.stats
        cpx_txt=${bam_file}.CpX.stats.txt
        mv ${bam_file}_stat.txt $cpx_txt

        ## convert BED to BigBED using bedToBigBed
        print_str "Converting BED ($cpg_bed) to BigBed"
        bed_bb=${cpg_bed%.bed}.bb
        tmp_file=$(/src/mcall_bed_reformat.pl $bed_in)
        bedToBigBed $tmp_file ${chrom_sizes} $bed_bb
        #rm -f $tmp_file

    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${preemtible}"
    }
    output {
        ???
    }
}

#https://multiqc.info/docs/
### Should have access to:
#fastqc for raw fastqs
#fastqc for trimmed fastqs
#cutadapt
task multiqc {
    File fastq1
    File fastq2
    File fastq1_trimmed
    File fastq2_trimmed
    File? bam_markdups_report

    command {
        mkdir multiqc_indir
        mv ${bam} multiqc_indir/
        mv ${fastq1} multiqc_indir/
        mv ${fastq2} multiqc_indir/

        multiqc multiqc_indir --filename ${sample_id}.multiqc_report.txt
    }

    runtime {
        docker: "${docker}"
        memory: "${mem}GBB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${preemtible}"
    }
    output {
        multiqc_report_html="${sample_id}.multiqc_report.txt"
    }
}

workflow bsmap_to_mcall_PE {

    String? docker="adunford/bsmap_to_mcall:0.32"

    call preprocess_fastqs {
    input:
        fastq1=fastq1
        fastq2=fastq2
        docker=docker
    }

    call fastqc as fastqc_raw {

    }

    call fastqc as fastqc_preprocessed {

    }

    call bsmap{

    }

    if(run_clip_overlap){
        call clip_overlaps
    }

    if(run_markduplicates){
        call markduplicates{
            input:
                bam_file = select_first([run_clip_overlap.bam_overlap_clipped , bsmap.bam])
                bam_index = select_first([run_clip_overlap.bam_overlap_clipped_index , bsmap.bam_index])
        }
    }

    call mcall {
        input:
            bam_file = select_first([markduplicates.bam_md, run_clip_overlap.bam_overlap_clipped , bsmap.bam])
            bam_index = select_first([markduplicates.bam_md_index, run_clip_overlap.bam_overlap_clipped_index , bsmap.bam_index])
    }

    call multiqc{
        ### Get .zip & .html for fastqc ; trimmed report; markdup report
        ### Outputs will include a tar.gz + multiqc output
        input:
            ??? #fastqc - raw
            ??? #fastqc - processed
            ??? #cutadapt log
            ??? #bsmap bam - raw
            ??? #bsmap bam - post-markduplicates ?
    }

    output {
	    File bsmap_bam_final = select_first([markduplicates.bam_md, clip_overlaps.bam_overlap_clipped, bsmap.bam])
        Array[File] bsmap_align_stats = select_all([bsmap.bam_stats, clip_overlaps.bam_overlap_clipped_stats, markduplicates.bam_md_stats])
        #File mcall ???
        #File multiqc ???
    }
}
