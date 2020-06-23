task bsmap_to_mcall{

	#FQ1       : first fastq file of a pair
	File fastq1
	#FQ2       : second fastq file of a pair
	File fastq2
	#SAMPLE    : sample name, used for output naming
	String sample_id
	#GENOME    : either 'mm9'/'mm10' (mouse) or 'hg19'/'hg38' (human).
	String genome_reference
	#C_THREADS : computing(alignment) threads to use.
	Int threads
	#n         : special flag for 2-color devices NextSeq and NovaSeq.
	#a         : perform not only CpG, but additionally CpA,C,T calling
	#o         : overlap-based clipping using bamutils (experimental)

  	Int disk_size
  	Int mem_size


  	command <<<
		sh /src/WGBS_Workflow_01_PE.sh
 	>>>

  	runtime {
    	docker: "adunford/bsmap_to_mcall:0.1"
    	memory: mem_size
    	disks: "local-disk " + disk_size + " HDD"
    	preemptible: 3
	}
  	output {
    	File MAF_Logs_tonly = "${case_name}.GermSomLogodds.maf"
    	File MAF_Pass_muts_tonly = "${case_name}.GermSomLogodds.pass.maf"
    	File Seg_file_Pass_tonly = "${case_name}.GermSomLogodds.seg.txt"
  	}
}

workflow tonly_pipeline {
  	call bsmap_to_mcall_PE
}