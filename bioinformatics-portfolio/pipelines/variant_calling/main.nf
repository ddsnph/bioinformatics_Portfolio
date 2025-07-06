#!/usr/bin/env nextflow

params.reads1     = "inputs/SRR396637.sra_1.fastq"
params.reads2     = "inputs/SRR396637.sra_2.fastq"
params.reference  = "inputs/pao1.fa"
params.outdir     = "results"
params.logdir     = "logs"

process INDEX_REF {
    tag "Index reference"
    input:
      file ref from params.reference
    output:
      file "*.bwt" into refidx
    script:
    """
    bwa index $ref > ${params.logdir}/bwa_index.log 2>&1
    """
}

process ALIGN {
    tag "Align reads"
    input:
      file ref from params.reference
      file r1 from params.reads1
      file r2 from params.reads2
    output:
      file "aligned.bam" into alignbam
    script:
    """
    bwa mem $ref $r1 $r2 | samtools view -Sb - > aligned.bam 2> ${params.logdir}/bwa_mem.log
    """
}

process SORT_BAM {
    tag "Sort BAM"
    input:
      file bam from alignbam
    output:
      file "sorted.bam" into sortedbam
    script:
    """
    samtools sort -o sorted.bam $bam 2> ${params.logdir}/samtools_sort.log
    """
}

process INDEX_BAM {
    tag "Index BAM"
    input:
      file bam from sortedbam
    output:
      file "sorted.bam.bai"
    script:
    """
    samtools index $bam 2> ${params.logdir}/samtools_index.log
    """
}

process CALL_VARIANTS {
    tag "Call variants"
    input:
      file bam from sortedbam
      file ref from params.reference
    output:
      file "variants.vcf"
    script:
    """
    bcftools mpileup -Ou -f $ref $bam 2> ${params.logdir}/bcftools_mpileup.log | \
    bcftools call -mv -Oz -o variants.vcf.gz 2> ${params.logdir}/bcftools_call.log
    bcftools index variants.vcf.gz
    bcftools view variants.vcf.gz > variants.vcf
    """
}

workflow {
    INDEX_REF()
    ALIGN()
    SORT_BAM()
    INDEX_BAM()
    CALL_VARIANTS()
}
