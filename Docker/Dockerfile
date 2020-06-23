from broadinstitute/gatk:4
#author adunford

#Install prerequisites
#pigz
RUN apt-get update -y
RUN apt-get install -y pigz
#gatk (Already installed)

#bsmap
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install bsmap


#mcall
RUN conda install -c bioconda moabs

#fastqc
RUN apt install -y default-jre
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
RUN unzip fastqc_v0.11.9.zip

#multiqc
RUN conda install -c bioconda -c conda-forge multiqc

#Set to POSIX in the script currently being used...
#ENV LC_ALL "C.UTF-8"
#ENV LANG "C.UTF-8" 
ENV LC_ALL "POSIX"
ENV LANG "POSIX"

#bamutils
RUN conda install bamutil

#cutadapt
RUN conda install -c bioconda cutadapt

#bamCoverage
RUN conda install -c bioconda deeptools




RUN mkdir /src
#need to update samtools to later version than the one gatk has
RUN conda install -c bioconda samtools

#update tools useful for development
RUN apt-get install -y emacs tmux

#custom scripts
COPY WGBS_Workflow_01_PE_AD.sh WGBS_Workflow_01_SR.sh bamStat.pl /src/
RUN chmod +x /src/*
