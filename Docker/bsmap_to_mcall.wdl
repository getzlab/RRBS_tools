task bsmap_to_mcall{

    #FQ1       : first fastq file of a pair
    File fastq1
    #FQ2       : second fastq file of a pair
    File fastq2
    #SAMPLE    : sample name, used for output naming
    String sample_id
    #GENOME    : either 'mm9'/'mm10' (mouse) or 'hg19'/'hg38' (human).
    String genome_reference
    File reference_fa
    File reference_sizes
    #C_THREADS : computing(alignment) threads to use.
    Int threads
    Int seed_size #seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.
    #n         : special flag for 2-color devices NextSeq and NovaSeq.
    #a         : perform not only CpG, but additionally CpA,C,T calling
    #o         : overlap-based clipping using bamutils (experimental)
    String docker_tag
    Int disk_size
    Int mem_size
    Int preemtible


    command <<<
    	    /src/WGBS_Workflow_01_PE_AD.sh -1 ${fastq1} -2 ${fastq2} -s ${sample_id} -g ${genome_reference} -r ${reference_fa} -z ${reference_sizes}  -c ${threads} -q ${seed_size}
 	    >>>

    runtime {
    	    docker: "adunford/bsmap_to_mcall:" + docker_tag
    	    memory: mem_size
            CPU: threads +1
    	    disks: "local-disk " + disk_size + " HDD"
    	    preemptible: preemtible
    }
    output {
	    File results_tar = "results.tar"
    }
}

workflow bsmap_to_mcall_PE {
  	 call bsmap_to_mcall
}