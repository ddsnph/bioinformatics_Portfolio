#!/usr/bin/env nextflow
/*
========================================================================================
 VARIANT CALLING PIPELINE
========================================================================================
 Author:      Myren DeShawn Sutton
 Date:        July 6, 2025
 Description: A Nextflow pipeline for identifying genetic variants from paired-end
              sequencing reads. The workflow performs reference indexing, read
              alignment, sorting, indexing, and variant calling.
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

params.reads     = "data/*_{1,2}.fastq"
params.reference = "data/pao1.fa"
params.outdir    = "results"
params.profile   = "standard"

params.s3_reads     = "s3://bioinformatics-demo-july2025/inputs/*_{1,2}.fastq"
params.s3_reference = "s3://bioinformatics-demo-july2025/inputs/pao1.fa"
params.s3_outdir    = "s3://bioinformatics-demo-july2025/results"

log.info """
         V A R I A N T - C A L L I N G - P I P E L I N E
         =============================================
         Reference: ${params.reference}
         Reads:     ${params.reads}
         Output:    ${params.outdir}
         """

workflow {
    def reads_path = (params.profile && params.profile.contains('aws')) ? params.s3_reads : params.reads
    def reference_path = (params.profile && params.profile.contains('aws')) ? params.s3_reference : params.reference
    def output_dir = (params.profile && params.profile.contains('aws')) ? params.s3_outdir : params.outdir
    
    ref_ch = Channel.fromPath(reference_path)
    read_pairs_ch = Channel.fromFilePairs(reads_path)

    INDEX_REF(ref_ch)
    ALIGN(read_pairs_ch, INDEX_REF.out)
    SORT_BAM(ALIGN.out)
    INDEX_BAM(SORT_BAM.out)
    CALL_VARIANTS(SORT_BAM.out, INDEX_BAM.out, ref_ch)
}

process INDEX_REF {
    queue 'nextflow-queue'
    
    input:
      path reference
    output:
      path "bwa_index", emit: bwa_index
    script:
    """
    mkdir bwa_index
    bwa index -p bwa_index/pao1 ${reference}
    """
}

process ALIGN {
    queue 'nextflow-queue'
    
    input:
      tuple val(sample_id), path(reads)
      path bwa_index_dir
    output:
      tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    script:
    """
    bwa mem ${bwa_index_dir}/pao1 ${reads[0]} ${reads[1]} | samtools view -Sb - > ${sample_id}.bam
    """
}

process SORT_BAM {
    queue 'nextflow-queue'
    
    input:
      tuple val(sample_id), path(bam_file)
    output:
      tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: sorted_bam
    script:
    """
    samtools sort ${bam_file} -o ${sample_id}.sorted.bam
    """
}

process INDEX_BAM {
    queue 'nextflow-queue'
    
    input:
      tuple val(sample_id), path(sorted_bam_file)
    output:
      tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: bai
    script:
    """
    samtools index ${sorted_bam_file}
    """
}

process CALL_VARIANTS {
    queue 'nextflow-queue'
    publishDir "${params.outdir}", mode: 'copy'

    input:
      tuple val(sample_id), path(sorted_bam)
      tuple val(sample_id2), path(bai)
      path reference

    output:
      tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf
      tuple val(sample_id), path("${sample_id}.vcf.gz.tbi"), emit: tbi

    script:
    """
    samtools faidx ${reference}
    bcftools mpileup -f ${reference} ${sorted_bam} | bcftools call -mv -Oz -o ${sample_id}.vcf.gz
    tabix -p vcf ${sample_id}.vcf.gz
    """
}

