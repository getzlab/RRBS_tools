from continuumio/miniconda3:4.8.2
#author adunford

#Install prerequisites
#pigz, random unix tools
RUN mkdir -p /usr/share/man/man1
RUN apt-get update -y
RUN apt-get install -y pigz tmux emacs curl wget  default-jre

##Download and install conda
#WORKDIR /tmp
#RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
#RUN  yes | bash Anaconda3-2020.02-Linux-x86_64.sh
#RUN source ~/.bashrc
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

#gatk
#RUN conda install -c bioconda gatk4

#bsmap, gatk5, fastqc, multiqc and other stuff I haven't had a reason to install separately with conda
RUN conda install bsmap gatk4 fastqc  multiqc bamutil cutadapt deeptools  ucsc-bedtobigbed 


#mcall
RUN conda install -c bioconda moabs

#fastqc
#RUN apt install -y default-jre
#RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
#RUN unzip fastqc_v0.11.9.zip

#multiqc
#RUN conda install -c bioconda -c conda-forge multiqc

#Set to POSIX in the script currently being used...
ENV LC_ALL "POSIX"
ENV LANG "POSIX"

#bamutils
#RUN conda install bamutil

#cutadapt
#RUN conda install -c bioconda cutadapt

#bamCoverage
#RUN conda install -c bioconda deeptools




RUN mkdir /src
#need to update samtools to later version for whatever reason
RUN conda install -c bioconda samtools


#custom scripts
COPY WGBS_Workflow_01_PE_AD.sh WGBS_Workflow_01_SR.sh *.pl /src/
RUN chmod +x /src/*

