task bismark_pe {
    File fastq1
    File fastq2
    File genome_index
    File chrom_sizes
    String sample_id
    Int multicore=1 #multicore=1 is equivalent to using 3-4 CPUs for parallelizable bismark tasks
    String? trim_args=''
    String? bismark_args=''

    Int? memory = "24"
    Int? num_threads = "4"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "2"
    Int? disk_size_gb = ceil( (size(fastq1, "G") + size(fastq2, "G")) * disk_scaler) + disk_size_buffer
    String docker = "gcr.io/broad-cga-bknisbac-wupo1/bismark:0.1"
    String? disk_type = "SSD"
    Int? num_preempt = "4"

    command {
        set -euo pipefail
        ### start monitoring script
        /src/monitor_script.sh > ${sample_id}.monitoring.log &

        ### trim fastqs
        trim_galore --paired ${trim_args} ${fastq1} ${fastq2}
        mv *.1_val_1.fq.gz ${sample_id}.1_val_1.fastq.gz
        mv *.2_val_2.fq.gz ${sample_id}.2_val_2.fastq.gz
        echo Printing dir contents before fq.gz del:
        find . -type f | xargs ls -l #print files names to output
        rm *.fq.gz #remove intermediate fastq to minimize required disk space

        ### bismark alignment
        mkdir bismark_index
        tar zxf ${genome_index} -C bismark_index

        bismark ${bismark_args} --genome bismark_index --multicore ${multicore} --basename ${sample_id} -1 ${sample_id}.1_val_1.fastq.gz -2 ${sample_id}.2_val_2.fastq.gz

        echo Printing dir contents after bismark_pe alignment:
        find . -type f | xargs ls -l #print files names to output

        #reformat file names for convenience
        mv *PE_report.txt ${sample_id}_report.txt
        mv *_pe.bam ${sample_id}.bam

        #sort and index bam
        samtools sort -o ${sample_id}.sorted.bam ${sample_id}.bam -@ ${num_threads}
        samtools index ${sample_id}.sorted.bam ${sample_id}.sorted.bai

        ### bismark report
        bismark2report --alignment_report ${sample_id}_report.txt --output ${sample_id}_bismark_report.html

        ### get methylation scores
        bismark_methylation_extractor --multicore ${multicore} --gzip --bedGraph --buffer_size 50% --genome_folder bismark_index ${sample_id}.bam
        #convert bedgraph to bw
        gunzip -c "${sample_id}.bedGraph.gz" | bedtools sort > ${sample_id}.sorted.bedGraph
        /usr/local/bin/bedGraphToBigWig ${sample_id}.sorted.bedGraph ${chrom_sizes} ${sample_id}.bw

        echo Printing dir contents at end of bismark_pe task:
        find . -type f | xargs ls -l #print file names to output
    }
    runtime {
        continueOnReturnCode: false
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_size_gb} ${disk_type}"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
    output {
        File output_bam = "${sample_id}.sorted.bam"
        File output_bai = "${sample_id}.sorted.bai"
        File output_covgz = "${sample_id}.bismark.cov.gz"
        File report = "${sample_id}_report.txt"
        File bismark_report_html = "${sample_id}_bismark_report.html"
        File monitoring_log = "${sample_id}.monitoring.log"
        File mbias_report = "${sample_id}.M-bias.txt"
        File output_bigwig = "${sample_id}.bw"
    }
}


workflow bismark_pe_workflow {
    call bismark_pe
}
