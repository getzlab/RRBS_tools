task bamstats {
    File bam_file
    File bam_index
    String prefix = basename(bam_file, ".bam")

    #runtime inputs
    String? docker="gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.1"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_buffer = "10"
    Int? disk_scaler = "1"
    Int? disk_size_gb = ceil( (size(bam_file, "G") + size(bam_index, "G")) * disk_scaler) + disk_size_buffer
    Int? num_preempt = "4"

    command {
            /src/monitor_script.sh > ${prefix}.monitoring.log &
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
        File monitoring_log = "${prefix}.monitoring.log"
    }
}

workflow bamstats_workflow {
    call bamstats
}
