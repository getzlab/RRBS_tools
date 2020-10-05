task fastqc{
    #inputs from workflow config, or upstream task if preprocessing done
    File fastq1
    File fastq2
    String? fastq_suffix = '.fq.gz'
    String? fastq1_prefix = basename(basename(basename(basename(fastq1, fastq_suffix), '.gz'), '.fastq'), '.fq')
    String? fastq2_prefix = basename(basename(basename(basename(fastq2, fastq_suffix), '.gz'), '.fastq'), '.fq')

    #String? tar_gz_prefix="fastqc"
    String? fastqc_args = ""

    #runtime inputs
    String? docker = "adunford/bsmap_to_mcall:0.32"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "1"
    Int? disk_size_gb = ceil( (size(fastq1, "G") + size(fastq2, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh &

        mkdir fastqc_results

        fastqc \
        --outdir .
          --threads ${threads} \
          ${fastqc_args} \
          ${fastq1} ${fastq2}

          find . | xargs ls -l
    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File fq1_zip="${fastq1_prefix}_fastqc.zip"
        File fq2_zip="${fastq2_prefix}_fastqc.zip"
        File fq1_html="${fastq1_prefix}_fastqc.html"
        File fq2_html="${fastq2_prefix}_fastqc.html"
    }
}

task trim_fastqs{
    #inputs from workflow config
    File fastq1
    File fastq2
    String sample_id #for log output

    String? fastq_suffix = ".fq.gz"
    #Don't need to expose these 3
    String? fastq1_prefix = basename(fastq1, fastq_suffix)
    String? fastq2_prefix = basename(fastq1, fastq_suffix)
    String trimmed_suffix = ".trm.fq.gz"
    # adaptors[illumina]='AGATCGGAAGAGC'
    # adaptors[smallrna]='TGGAATTCTCGG'
    # adaptors[nextera]='CTGTCTCTTATA'
    String? adaptors="AGATCGGAAGAGC" #illumina
    Int? min_length="25"
    Int? min_qual="20"
    Int? swift_trim_off="10"
    # If data has been generated on NextSeq or NovaSeq 2-color machines,
    # we need to apply a different type of quality clipping (due to dark cycles)
    # If so, specify "--nextseq-trim"
    String? trimming_type="--quality-cutoff"

    #runtime inputs
    String? docker="adunford/bsmap_to_mcall:0.32"
    Int? mem = "3"
    Int? threads = "2"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "2"
    Int? disk_size_gb = ceil( (size(fastq1, "G") + size(fastq2, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh &

            ### cutadapt
            ## R1----->
            ## ~~~~~~~|=======================|~~~~~~~
            ##   AD   |        DNA            |   AD
            ## ~~~~~~~|=======================|~~~~~~~
            ##                                <-----R2
            echo current adaptor sequences: ${adaptors}
            # quality clipping is done BEFORE adaptor clipping
            cutadapt \
        	${trimming_type} ${min_qual} \
        	--overlap 5 \
        	--minimum-length $((${min_length} + 2 * ${swift_trim_off})) \
        	--cores $((${threads} / 2)) \
        	--adapter ${adaptors} \
        	-A ${adaptors} \
        	--interleaved \
        	${fastq1} \
        	${fastq2} \
        	2>.log_a | \
            cutadapt \
                --interleaved \
        	--minimum-length ${min_length} \
        	--cores $((${threads} / 2)) \
        	--cut ${swift_trim_off} \
        	--cut -${swift_trim_off} \
        	-U ${swift_trim_off} -U -${swift_trim_off} \
        	--output ${fastq1_prefix}${trimmed_suffix} \
        	--paired-output ${fastq2_prefix}${trimmed_suffix} \
        	- \
            &>.log_b

            ## merge log output
            cat .log_a .log_b > ${sample_id}_FastQ-Trimming-Adp.log

            find . | xargs ls -l
 	    }

    runtime {
    	    docker: "${docker}"
    	    memory: "${mem}GB"
            cpu: "${threads}"
    	    disks: "local-disk ${disk_size_gb} HDD"
    	    preemptible: "${num_preempt}"
    }
    output {
	    File fastq1_trimmed = "${fastq1_prefix}${trimmed_suffix}"
	    File fastq2_trimmed = "${fastq2_prefix}${trimmed_suffix}"
        File trimming_log = "${sample_id}_FastQ-Trimming-Adp.log"
    }
}

task bamstat {
    File bam_file
    File bam_index
    String prefix = basename(bam_file, ".bam")

    #runtime inputs
    String? docker="adunford/bsmap_to_mcall:0.32"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "1"
    Int? disk_size_gb = ceil( (size(bam_file, "G") + size(bam_index, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
            /src/monitor_script.sh &
            bamstats --bam ${bam_file} --threads $((${threads}-1)) --verbose | tee ${prefix}.stats.txt

            find . | xargs ls -l
        }

    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File bam_stats = "${prefix}.stats.txt"
    }
}

task bsmap{
    File fastq1
    File fastq2
    String sample_id
    File reference_fa
    Int? seed_size="12" #default=12(RRBS mode), 16(WGBS mode). min=8, max=16.
    Int? max_insert_size="1000" # max insert size for PE mapping (-x)
    # -q is quality threshold.  Here it's 20, default is 0, should we do 0 if preprocessing done?
    # -w<int>   maximum number of equal best hits to count, <=1000
    # -S seed for rng.  0 for system clock (not reproducible) otherwise produces reproducible results.
    # -u   report unmapped reads, default=off
    # -R          print corresponding reference sequences in SAM output, default=off
    String? bsmap_args="-q 20 -w 100 -S 1 -u -R" #??? understand each param in context of RRBS vs WGBS

    Int? sort_threads = "2"
    String? sort_args = ""

    #runtime inputs
    String? docker="adunford/bsmap_to_mcall:0.32"
    #mem=10, threads=12 is sufficient for "samtools sort" default mem/thread usage [768 MiB]
    Int? mem = "10"
    Int? threads = "12"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "2"
    Int? disk_size_gb = ceil( (size(fastq1, "G") + size(fastq2, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
            /src/monitor_script.sh &

            # -s = seed size, default=12(RRBS mode), 16(WGBS mode). min=8, max=16
            # -u = report unmapped reads, default=off
            # -R = print corresponding reference sequences in SAM output, default=off
            # -q =  quality threshold in trimming, 0-40, default=0 (no trim)
            # -w =  maximum number of equal best hits to count, <=1000
            # -S = 1 seed for random number generation used in selecting multiple hits, keep nonzero for reproucibility
            bam_file_unsorted=${sample_id}.bsmap.unsorted.bam
            bam_file=${sample_id}.bsmap.srt.bam

            bsmap \
              -v 0.1 -s ${seed_size} ${bsmap_args} \
              -x ${max_insert_size} \
              -p ${threads} \
              -d ${reference_fa} \
              -a ${fastq1} \
              -b ${fastq2} > $bam_file_unsorted

            samtools sort ${sort_args} --threads ${sort_threads} --output-fmt BAM -o $bam_file $bam_file_unsorted
            rm $bam_file_unsorted

            echo "Creating BAM Index for $bam_file"
            samtools index -@ ${threads} $bam_file

            find . | xargs ls -l
 	    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File bam = "${sample_id}.bsmap.srt.bam"
        File bam_index = "${sample_id}.bsmap.srt.bam.bai"
    }
}

task clip_overlap {
    #From upstream tasks
    File bam_file
    File bam_index
    String? bam_prefix=basename(bam_file, ".bam")
    String? clip_overlap_args="--stats --params --unmapped --noPhoneHome --poolSize 10000000"

    #runtime inputs
    String? docker="adunford/bsmap_to_mcall:0.32"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "20"
    Int? disk_scaler = "2"
    Int? disk_size_gb = ceil( size(bam_file, "G") + size(bam_index, "G") * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh &

        $bam_file_oc=${bam_prefix}.oc.bam

        bamutils clipOverlap --in ${bam_file} --out $bam_file_oc ${clip_overlap_args}

        echo $(date +"### [%Y-%m-%d %H:%M:%S] Creating BAM Index for $bam_file_oc")
        samtools index -@ ${threads} $bam_file_oc

        find . | xargs ls -l
    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }

    output {
        File bam_overlap_clipped="${bam_prefix}.oc.bam"
        File bam_overlap_clipped_index="${bam_prefix}.oc.bam.bai"
    }
}

task markduplicates {
    File bam_file
    File bam_index
    String sample_id

    String? bam_prefix=basename(bam_file, ".bam")

    String? remove_dups="false" #WARNING: "true" LEADS TO LOSS OF READS IN FINAL BAM
    Int? max_records_in_ram="500000" #As in GTEx; Helene: "10000000"
    Int? max_open_files="8000" #GATK default
    Int? java_mem="3" # max heap size #As in GTEx; Helene: "16G"
    String? read_name_regex="null" #skips optical duplicate finding (this is setting for HiSeq 4000 and NovaSeq 6000 that shouldnt have optical duplicates)
    String? other_args=""

    #runtime inputs
    String? docker = "adunford/bsmap_to_mcall:0.32"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "20"
    Int? disk_scaler = "3" # 3 for potential intermediate files
    Int? disk_size_gb = ceil( size(bam_file, "G") + size(bam_index, "G") * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh &

        bam_dedup=${bam_prefix}.md.bam

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
          --OUTPUT=${bam_prefix}.md.bam \
          --METRICS_FILE=${bam_prefix}.md-metrics.txt \ #if name not essential better to change to md-metrics.txt
          --TMP_DIR=$temp_dir \
          --VALIDATION_STRINGENCY=LENIENT \
          --READ_NAME_REGEX=${read_name_regex} \
          ${other_args}

        echo $(date +"### [%Y-%m-%d %H:%M:%S] Creating BAM Index for ${bam_prefix}.md.bam")
        samtools index -@ ${threads} ${bam_prefix}.md.bam

        find . | xargs ls -l
    }

    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File bam_md="${bam_prefix}.md.bam"
        File bam_md_index="${bam_prefix}.md.bam.bai"
        File bam_md_metrics="${bam_prefix}.md-metrics.txt"
    }
}

task mcall {
    File bam_file #either bsmap_bam or bam_md
    File bam_index
    String sample_id
    File reference_fa
    File reference_sizes

    String? mcall_args=""

    String? bam_prefix=basename(bam_file, ".bam")

    #runtime inputs
    String? docker = "adunford/bsmap_to_mcall:0.32"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "2"
    Int? disk_size_gb = ceil( size(bam_file, "G") + size(bam_index, "G") * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh &

        echo $(date +"### [%Y-%m-%d %H:%M:%S] Running Methylation Calling CpG")

        mcall \
          --threads ${threads} \
          --reference ${reference_fa} \
          --sampleName ${sample_id} \ #BK: removed genome index
          --mappedFiles ${bam_file} \
          ${mcall_args}

        ## output generated by mcall
        ## run.config.64964
        ## AM-RRBS-88_mm9.G.bed -> AM-WGBS-88_mm9.bsmap.srt.oc.rd.bam.G.bed
        ## AM-WGBS-88_mm9.bsmap.srt.oc.rd.bam.G.bed
        ## AM-WGBS-88_mm9.bsmap.srt.oc.rd.bam.HG.bed
        ## AM-WGBS-88_mm9.bsmap.srt.oc.rd.bam_stat.txt
        ## ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        ## <---        $bam_file         --->
        mv ${bam_file}.G.bed ${sample_id}.CpG.bed

        #Stats are not unique to CpG (would be the same for e.g. CpA analysis), hence naming CpX.stats
        mv ${bam_file}_stat.txt ${sample_id}.CpX.stats.txt

        ## convert BED to BigBED using bedToBigBed
        echo "Converting BED (${sample_id}.CpG.bed) to BigBed"
        tmp_file=$(/src/mcall_bed_reformat.pl ${sample_id}.CpG.bed)
        bedToBigBed $tmp_file ${reference_sizes} ${sample_id}.CpG.bb
        rm -f $tmp_file

        gzip ${sample_id}.CpG.bed

        find . | xargs ls -l
    }
    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File mcall_cpg_bed_gz = "${sample_id}.CpG.bed.gz"
        File mcall_cpg_bigbed = "${sample_id}.CpG.bb"
        File mcall_stats = "${sample_id}.CpX.stats.txt"
    }
}

#https://multiqc.info/docs/
### Should have access to:
#fastqc for raw fastqs
#fastqc for trimmed fastqs
#cutadapt
task multiqc {
    String sample_id
    File fastq1_fastqc
    File fastq2_fastqc
    File? fastq1_trimmed_fastqc
    File? fastq2_trimmed_fastqc
    File? trimming_log
    File? bam_markdups_report
    String? multiqc_args = ""

    #runtime inputs
    String? docker = "adunford/bsmap_to_mcall:0.32"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "1"
    Int? disk_size_gb = ceil( (size(fastq1_fastqc, "G") + size(fastq2_fastqc, "G") + size(fastq1_trimmed_fastqc, "G") + size(fastq2_trimmed_fastqc, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
        /src/monitor_script.sh &

        $multiqc_dir="${sample_id}.multiqc"
        mkdir $multiqc_dir

        mv ${fastq1_fastqc} $multiqc_dir/
        mv ${fastq2_fastqc} $multiqc_dir/

        if [ -n "${trimming_log}" ]; then
            mv ${trimming_log} $multiqc_dir/
        fi

        if [ -n "${fastq1_trimmed_fastqc}" ]; then
            mv ${fastq1_trimmed_fastqc} $multiqc_dir/
        fi

        if [ -n "${fastq2_trimmed_fastqc}" ]; then
            mv ${fastq2_trimmed_fastqc} $multiqc_dir/
        fi

        if [ -n "${bam_markdups_report}" ]; then
            mv ${bam_markdups_report} $multiqc_dir/
        fi

        multiqc $multiqc_dir --filename ${sample_id}.multiqc_report.txt ${multiqc_args}

        tar -cxvf $multiqc_dir.tar.gz $multiqc_dir

        find . | xargs ls -l
    }

    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File multiqc_report_html="${sample_id}.multiqc_report.txt"
        File multiqc_tar_gz="${sample_id}.multiqc.tar.gz"
    }
}

workflow bsmap_to_mcall_PE {
    String sample_id
    File fastq1
    File fastq2
    File reference_fa #for bsmap, mcall
    File reference_sizes #for bsmap, mcall

    ## Control variables
    Boolean? run_clip_overlap = true
    Boolean? run_markduplicates = true

    ### workflow-wide optional runtime variables
    String? docker="adunford/bsmap_to_mcall:0.32"
    Int? num_preempt="4"

    #### Per tasks ####
    ## for fastqc
    String? fastq_suffix #defaults to "fq.gz"
    String? fastqc_args_raw
    String? fastqc_args_trimmed

    ## for trim_fastqs
    String? trim_fastqs_fastq_suffix
    String? trim_fastqs_adaptors
    Int? trim_fastqs_min_length
    Int? trim_fastqs_min_qual
    Int? trim_fastqs_swift_trim_off
    String? trim_fastqs_trimming_type

    ## for bsmap
    Int? bsmap_seed_size #default=12(RRBS mode), 16(WGBS mode). min=8, max=16.
    Int? bsmap_max_insert_size # max insert size for PE mapping (-x)
    String? bsmap_args
    Int? bsmap_sort_threads
    Int? bsmap_sort_args

    ## for markduplicates
    String? markdups_remove_dups #WARNING: "true" LEADS TO LOSS OF READS IN FINAL BAM
    Int? markdups_max_records_in_ram #As in GTEx; Helene: "10000000"
    Int? markdups_max_open_files #GATK default
    Int? markdups_java_mem # max heap size #As in GTEx; Helene: "16G"
    String? markdups_read_name_regex #skips optical duplicate finding (this is setting for HiSeq 4000 and NovaSeq 6000 that shouldnt have optical duplicates)
    String? markdups_other_args

    ## for clip overlaps
    String? clip_overlap_args

    ## for mcall
    String? mcall_args

    ## for multiqc
    String? multiqc_args

    ### Per-task optional runtime variables
    Int? trim_fastqs_mem
    Int? trim_fastqs_threads
    Int? trim_fastqs_disk_size_buffer
    Int? trim_fastqs_num_preempt

    Int? fastqc_mem
    Int? fastqc_threads
    Int? fastqc_disk_size_buffer
    Int? fastqc_num_preempt

    Int? bsmap_mem
    Int? bsmap_threads
    Int? bsmap_disk_size_buffer
    Int? bsmap_num_preempt

    Int? bamstat_mem
    Int? bamstat_threads
    Int? bamstat_disk_size_buffer
    Int? bamstat_num_preempt

    Int? clip_overlap_mem
    Int? clip_overlap_threads
    Int? clip_overlap_disk_size_buffer
    Int? clip_overlap_num_preempt

    Int? markdups_mem
    Int? markdups_threads
    Int? markdups_disk_size_buffer
    Int? markdups_num_preempt

    Int? mcall_mem
    Int? mcall_threads
    Int? mcall_disk_size_buffer
    Int? mcall_num_preempt

    Int? multiqc_mem
    Int? multiqc_threads
    Int? multiqc_disk_size_buffer
    Int? multiqc_num_preempt


    call trim_fastqs {
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            sample_id=sample_id,
            docker = docker,
            mem = trim_fastqs_mem,
            threads = trim_fastqs_threads,
            disk_size_buffer = trim_fastqs_disk_size_buffer,
            num_preempt = select_first([trim_fastqs_num_preempt, num_preempt])
    }

    call fastqc as fastqc_raw {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            fastq_suffix = fastq_suffix,
            fastqc_args = fastqc_args_raw,
            docker = docker,
            mem = fastqc_mem,
            threads = fastqc_threads,
            disk_size_buffer = fastqc_disk_size_buffer,
            num_preempt = select_first([fastqc_num_preempt, num_preempt])
    }

    call fastqc as fastqc_trimmed {
        input:
            fastq1 = trim_fastqs.fastq1_trimmed,
            fastq2 = trim_fastqs.fastq2_trimmed,
            docker = docker,
            mem = fastqc_mem,
            threads = fastqc_threads,
            disk_size_buffer = fastqc_disk_size_buffer,
            num_preempt = select_first([fastqc_num_preempt, num_preempt])
    }

    File? trimmed_fastq1 = trim_fastqs.fastq1_trimmed
    File? trimmed_fastq2 = trim_fastqs.fastq2_trimmed

    call bsmap{
        input:
            fastq1=select_first([trimmed_fastq1, fastq1]),
            fastq2=select_first([trimmed_fastq2, fastq2]),
            sample_id=sample_id,
            reference_fa=reference_fa,
            seed_size=bsmap_seed_size,
            max_insert_size=bsmap_max_insert_size,
            bsmap_args=bsmap_args,
            sort_threads=bsmap_sort_threads,
            sort_args=bsmap_sort_args,
            docker = docker,
            mem = bsmap_mem,
            threads = bsmap_threads,
            disk_size_buffer = bsmap_disk_size_buffer,
            num_preempt = select_first([bsmap_num_preempt, num_preempt])
    }

    call bamstat as bamstat_bsmap{
        input:
            bam_file = bsmap.bam,
            bam_index = bsmap.bam_index,
            docker = docker,
            mem = bamstat_mem,
            threads = bamstat_threads,
            disk_size_buffer = bamstat_disk_size_buffer,
            num_preempt = select_first([bamstat_num_preempt, num_preempt])
    }

    if(run_clip_overlap==true){
        call clip_overlap{
            input:
                bam_file = bsmap.bam,
                bam_index = bsmap.bam_index,
                clip_overlap_args = clip_overlap_args,
                docker = docker,
                mem = clip_overlap_mem,
                threads = clip_overlap_threads,
                disk_size_buffer = clip_overlap_disk_size_buffer,
                num_preempt = select_first([clip_overlap_num_preempt, num_preempt])
        }
        call bamstat as bamstat_co{
            input:
                bam_file = clip_overlap.bam_overlap_clipped,
                bam_index = clip_overlap.bam_overlap_clipped_index,
                docker = docker,
                mem = bamstat_mem,
                threads = bamstat_threads,
                disk_size_buffer = bamstat_disk_size_buffer,
                num_preempt = select_first([bamstat_num_preempt, num_preempt])
        }
    }

File? overlap_clipped_bam = clip_overlap.bam_overlap_clipped
File? overlap_clipped_bam_index = clip_overlap.bam_overlap_clipped_index

    if(run_markduplicates==true){
        call markduplicates{
            input:
                bam_file = select_first([overlap_clipped_bam , bsmap.bam]),
                bam_index = select_first([overlap_clipped_bam_index , bsmap.bam_index]),
                sample_id = sample_id,
                remove_dups = markdups_remove_dups,
                max_records_in_ram = markdups_max_records_in_ram,
                max_open_files = markdups_max_open_files,
                java_mem = markdups_java_mem,
                read_name_regex = markdups_read_name_regex,
                other_args = markdups_other_args,
                docker = docker,
                mem = markdups_mem,
                threads = markdups_threads,
                disk_size_buffer = markdups_disk_size_buffer,
                num_preempt = select_first([markdups_num_preempt, num_preempt])
        }
        call bamstat as bamstat_md{
            input:
                bam_file = markduplicates.bam_md,
                bam_index = markduplicates.bam_md_index,
                docker = docker,
                mem = bamstat_mem,
                threads = bamstat_threads,
                disk_size_buffer = bamstat_disk_size_buffer,
                num_preempt = select_first([bamstat_num_preempt, num_preempt])
        }
    }

File? markdup_bam = clip_overlap.bam_overlap_clipped
File? markdup_bam_index = clip_overlap.bam_overlap_clipped_index


    call mcall {
        input:
            bam_file = select_first([markdup_bam, overlap_clipped_bam , bsmap.bam]),
            bam_index = select_first([markdup_bam_index, overlap_clipped_bam_index , bsmap.bam_index]),
            sample_id = sample_id,
            reference_fa=reference_fa,
            reference_sizes=reference_sizes,
            mcall_args = mcall_args,
            docker = docker,
            mem = mcall_mem,
            threads = mcall_threads,
            disk_size_buffer = mcall_disk_size_buffer,
            num_preempt = select_first([mcall_num_preempt, num_preempt])
    }

    call multiqc{
        ### Get .zip & .html for fastqc ; trimmed report; markdup report
        ### Outputs will include a tar.gz + multiqc output
        input:
            sample_id = sample_id,
            fastq1_fastqc = fastqc_raw.fq1_zip,
            fastq2_fastqc = fastqc_raw.fq2_zip,
            fastq1_trimmed_fastqc = trim_fastqs.fastq1_trimmed,
            fastq2_trimmed_fastqc = trim_fastqs.fastq2_trimmed,
            trimming_log = trim_fastqs.trimming_log,
            bam_markdups_report = markduplicates.bam_md_metrics,
            multiqc_args = multiqc_args,
            docker = docker,
            mem = multiqc_mem,
            threads = multiqc_threads,
            disk_size_buffer = multiqc_disk_size_buffer,
            num_preempt = select_first([multiqc_num_preempt, num_preempt])
    }


    File? bamstats_bsmap_file = bamstat_bsmap.bam_stats
    File? bamstats_co_file = bamstat_co.bam_stats
    File? bamstats_md_file = bamstat_md.bam_stats

    meta {
        author: "Binyamin A. Knisbacher"
        email: "bknisbac@broadinstitute.org"
        author: "Andrew Dunford"
        email: "adunford@broadinstitute.org"
    }

    output {
        #bsmap
        File bsmap_bam_final = select_first([markduplicates.bam_md, clip_overlap.bam_overlap_clipped, bsmap.bam])
        Array[File] bsmap_align_stats = select_all([bamstats_bsmap_file, bamstats_co_file, bamstats_md_file])
        #multiqc
        File multiqc_report_html = multiqc.multiqc_report_html
        File multiqc_tar_gz = multiqc.multiqc_tar_gz
        #mcall
        File mcall_cpg_bed_gz = mcall.mcall_cpg_bed_gz
        File mcall_cpg_bigbed = mcall.mcall_cpg_bigbed
        File mcall_stats = mcall.mcall_stats
    }
}
