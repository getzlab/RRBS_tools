task parse_bismark_report {
    File bismark_report_html
    String sample_id
    String paired_str = "PE"

    Float? memory = "3.75"
    Int? num_threads = "1"
    Int? disk_size_gb = "50"
    String docker = "gcr.io/broad-cga-bknisbac-wupo1/bismark:0.2"
    String? disk_type = "HDD"
    Int? num_preempt = "4"

    command {
        set -euo pipefail
        ### start monitoring script
        /src/monitor_script.sh > ${sample_id}.monitoring.log &

        mv ${bismark_report_html} ${sample_id}_bismark_report.html

        ### Parse bismark report
        /src/extract_bismark_reads_reported.py ${sample_id}_bismark_report.html ${paired_str}

    }
    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_size_gb} ${disk_type}"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }
    output {
        File bismark_report_tsv = "${sample_id}_bismark_bamstats.tsv"
    }
}

workflow parse_bismark_report_workflow {
    call parse_bismark_report
}
