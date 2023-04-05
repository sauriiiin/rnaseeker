FROM rocker/r-ver:4.1.1

RUN apt-get update && \
    apt-get install -y curl wget unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages(c("tidyverse", "tximport", "DESeq2", "reactome.db", "org.Sc.sgd.db", "clusterProfiler"))'

RUN wget -O salmon.tar.gz https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz && \
    tar -zxvf salmon.tar.gz && \
    mv salmon-1.5.2_linux_x86_64 /opt/ && \
    rm salmon.tar.gz

WORKDIR /app

COPY nextflow.config .
COPY main.nf .
COPY generate_reads.R .
COPY diff_exp.R .
COPY go_kegg.R .

ENTRYPOINT ["nextflow", "run", "main.nf", "-params-file", "nextflow.config"]
