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

// --- INPUT/OUTPUT PARAMETERS --- //
params.reads     = "data/*_{1,2}.fastq" // Path pattern for paired-end sequencing reads
params.reference = "data/pao1.fa"         // Path to the reference genome FASTA file
params.outdir    = "results"                // Directory where all final results will be saved

// --- PIPELINE LOGO --- //
log.info """
         V A R I A N T - C A L L I N G - P I P E L I N E
         =============================================
         Reference: ${params.reference}
         Reads:     ${params.reads}
         Output:    ${params.outdir}
         """

//======================================================================================
//                                 WORKFLOW DEFINITION
//======================================================================================
workflow {
    ref_ch = Channel.fromPath(params.reference)
    read_pairs_ch = Channel.fromFilePairs(params.reads)

    INDEX_REF(ref_ch)
    ALIGN(read_pairs_ch, INDEX_REF.out)
    SORT_BAM(ALIGN.out)
    INDEX_BAM(SORT_BAM.out)
    CALL_VARIANTS(SORT_BAM.out, INDEX_BAM.out, ref_ch)
}

//======================================================================================
//                                     PROCESSES
//======================================================================================

process INDEX_REF {
    input:
      path reference
    output:
      path "bwa_index"
    script:
    """
    mkdir bwa_index
    bwa index -p bwa_index/pao1 ${reference}
    """
}

process ALIGN {
    input:
      tuple val(sample_id), path(reads)
      path bwa_index_dir
    output:
      path "${sample_id}.bam"
    script:
    """
    bwa mem ${bwa_index_dir}/pao1 ${reads[0]} ${reads[1]} | samtools view -Sb - > ${sample_id}.bam
    """
}

process SORT_BAM {
    input:
      path bam_file
    output:
      path "*.sorted.bam"
    script:
    """
    samtools sort ${bam_file} -o ${bam_file}.sorted.bam
    """
}

process INDEX_BAM {
    input:
      path sorted_bam_file
    output:
      path "*.bai"
    script:
    """
    samtools index ${sorted_bam_file}
    """
}

process CALL_VARIANTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
      path sorted_bam
      path bai
      path reference

    output:
      path "*.vcf.gz"
      path "*.vcf.gz.tbi", optional: true

    script:
    def sample_id = sorted_bam.baseName.toString() - '.sorted'
    """
    samtools faidx ${reference}
    bcftools mpileup -f ${reference} ${sorted_bam} | bcftools call -mv -Oz -o ${sample_id}.vcf.gz
    
    # Create tabix index for the VCF file
    tabix -p vcf ${sample_id}.vcf.gz
    """
}

