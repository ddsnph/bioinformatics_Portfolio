#!/usr/bin/env nextflow

params.reads = "$baseDir/inputs/*.fastq.gz"
params.reference = "$baseDir/inputs/chr21.fa"

process alignReads {
    input:
    file reads from params.reads
    file ref from params.reference

    output:
    file "aligned.bam" into aligned_bam

    """
    bwa index $ref
    bwa mem $ref $reads | samtools view -Sb - > aligned.bam
    """
}

