# Base image
FROM rocker/r-ver:4.1.1

# Install system dependencies
RUN apt-get update && \
    apt-get install -y wget curl git build-essential libcurl4-gnutls-dev libssl-dev libxml2-dev libz-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Java
RUN apt-get update && \
    apt-get install -y openjdk-11-jre && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz && \
    tar -xzf salmon-1.5.2_linux_x86_64.tar.gz && \
    mv salmon-1.5.2_linux_x86_64/bin/* /usr/local/bin/ && \
    rm -rf salmon-1.5.2_linux_x86_64.tar.gz salmon-1.5.2_linux_x86_64

# Install TrimGalore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.7.tar.gz -o trim_galore.tar.gz && \
    tar -xzf trim_galore.tar.gz && \
    mv TrimGalore-0.6.7/trim_galore /usr/local/bin/ && \
    rm -rf TrimGalore-0.6.7 trim_galore.tar.gz

# Install R and packages
RUN R -e 'install.packages(c("tidyverse", "tximport", "DESeq2", "reactome.db", "org.Sc.sgd.db", "clusterProfiler"))'

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
