# Base image
FROM rocker/r-ver:4.1.2
# FROM r-base:latest

# Install system dependencies
RUN apt-get update && \
    apt-get install -y wget curl git build-essential libcurl4-gnutls-dev libssl-dev libxml2-dev libz-dev libudunits2-dev cutadapt openjdk-11-jre gawk bison libc6 libglpk-dev libgomp1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
    
# Install R and packages
RUN R -e 'install.packages(c("tidyverse","BiocManager"))'
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("tximport", \
    "DESeq2", "reactome.db", "org.Sc.sgd.db", "enrichplot", "clusterProfiler", \
    "GO.db", "apeglm", ask = F))'

# Install Salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v0.12.0/salmon-0.12.0_linux_x86_64.tar.gz && \
    tar -xzf salmon-0.12.0_linux_x86_64.tar.gz && \
    mv salmon-0.12.0_linux_x86_64/bin/* /usr/local/bin/ && \
    rm salmon-0.12.0_linux_x86_64/lib/libm.so.6 && \
    rm salmon-0.12.0_linux_x86_64/lib/libgomp.so.1 && \
    mv salmon-0.12.0_linux_x86_64/lib/* /usr/local/lib/ && \
    rm -rf salmon-0.12.0_linux_x86_64.tar.gz salmon-0.12.0_linux_x86_64
    
# Install GLIB 2.29
#RUN wget -c https://ftp.gnu.org/gnu/glibc/glibc-2.29.tar.gz && \
#    tar -zxvf glibc-2.29.tar.gz && \
#    cd glibc-2.29 && \
#    mkdir build && \
#    cd build && \
#    ../configure --prefix=/opt/glibc && \
#    make && \
#    make install && \
#    rm -rf glibc-2.29.tar.gz

# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    cd FastQC/ && \
    chmod +x fastqc && \
    ln -s fastqc /usr/local/bin/fastqc

# Install TrimGalore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.7.tar.gz -o trim_galore.tar.gz && \
    tar -xzf trim_galore.tar.gz && \
    mv TrimGalore-0.6.7/trim_galore /usr/local/bin/ && \
    rm -rf TrimGalore-0.6.7 trim_galore.tar.gz

# Install Nextflow
RUN curl -fsSL https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

RUN mkdir /workdir
WORKDIR /workdir

# Copy pipeline files
COPY nextflow.config .
COPY main.nf .
COPY generate_reads.R .
COPY diff_exp.R .
COPY go_kegg.R .

# Run the pipeline
CMD ["nextflow", "run", "-c", "nextflow.config"]

ENTRYPOINT ["nextflow", "run", "main.nf"]