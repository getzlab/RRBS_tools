from continuumio/miniconda3:4.8.2
#author Binyamin Knisbacher (adapted from Andrew Dunford)

#Install prerequisites: pigz
RUN apt-get update -y
RUN apt-get install -y apt-utils curl wget pigz python3 python3-pip \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

##Download and install conda
#WORKDIR /tmp
#RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
#RUN  yes | bash Anaconda3-2020.02-Linux-x86_64.sh
#RUN source ~/.bashrc
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install bismark samtools trim-galore

#custom scripts
RUN mkdir /src
COPY monitor_script.sh /src/
RUN chmod +x /src/*
