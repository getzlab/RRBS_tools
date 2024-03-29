task bamstat_to_tsv {
    File bamstat_file
    String sample_id
    String? seq_type = "PE"
    String? bamstat_args = ""

    #runtime inputs
    String? docker="gcr.io/broad-cga-bknisbac-wupo1/bisulfite_tools:0.4"
    Int? mem = "3"
    Int? threads = "1"
    Int? disk_size_gb = "50"
    Int? num_preempt = "4"

    command {
            /src/parse_bamstats.py -f ${bamstat_file} -s ${sample_id} -o ${sample_id}.bamstat.tsv -t ${seq_type} ${bamstat_args}
        }

    runtime {
        docker: "${docker}"
        memory: "${mem}GB"
        cpu: "${threads}"
        disks: "local-disk ${disk_size_gb} HDD"
        preemptible: "${num_preempt}"
    }
    output {
        File bam_stats_tsv = "${sample_id}.bamstat.tsv"
    }
}

workflow bamstats_workflow {
    call bamstat_to_tsv
}
