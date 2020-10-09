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
#Specifying gatk 4.1.5 because sometime before 4.1.8 the syntax changes
RUN conda install bsmap gatk4=4.1.5.0 fastqc multiqc bamutil cutadapt deeptools ucsc-bedtobigbed
RUN conda install pandas

#mcall
RUN conda install -c anaconda boost
#RUN conda install  moabs
#RUN wget https://s3.amazonaws.com/deqiangsun/software/moabs-v1.3.0.src.x86_64_Linux.data.tar.gz
RUN mkdir /src
WORKDIR /src
RUN git clone https://github.com/sunnyisgalaxy/moabs.git
#RUN export PATH=PATH:/src/moabs/bin
#WORKDIR  moabs
#RUN apt-get install make
#RUN make
#RUN make install

#Set to POSIX in the script currently being used...
ENV LC_ALL "POSIX"
ENV LANG "POSIX"

#install needed perl modules
RUN conda install -c bioconda perl-statistics-descriptive


#need to update samtools to later version and force install openssl 1.0 to keep libraries neat
RUN conda install -c bioconda samtools openssl=1.0
ENV PATH="${PATH}:/src/moabs/bin"

#Need to symlink mysql to older version expected by bedtobigbed... seems to work so far, may need to install older version manually as above with openssl
RUN ln -s /opt/conda/lib/libmysqlclient.so.20  /opt/conda/lib/libmysqlclient.so.18

#custom scripts
COPY *.sh *.pl *.py /src/
RUN chmod +x /src/*

RUN mkdir /workdir
WORKDIR /workdir
