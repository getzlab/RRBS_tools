from continuumio/miniconda3:4.8.2
#author Binyamin Knisbacher (adapted from Andrew Dunford)

#Install prerequisites: pigz
RUN apt-get update -y
RUN apt-get install -y pigz

##Download and install conda
#WORKDIR /tmp
#RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
#RUN  yes | bash Anaconda3-2020.02-Linux-x86_64.sh
#RUN source ~/.bashrc
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install bismark samtools trim_galore

#Set to POSIX in the script currently being used...
# ENV LC_ALL "POSIX"
# ENV LANG "POSIX"

#need to update samtools to later version and force install openssl 1.0 to keep libraries neat
RUN conda install -c bioconda samtools openssl=1.0
ENV PATH="${PATH}:/src/moabs/bin"

#custom scripts
RUN mkdir /src
COPY ../Docker/monitor_script.sh /src/
RUN chmod +x /src/*

RUN mkdir /workdir
WORKDIR /workdir
