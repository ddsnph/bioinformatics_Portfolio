FROM ubuntu:22.04

# Install bioinformatics tools
RUN apt-get update && apt-get install -y \
    bwa \
    samtools \
    bcftools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data