#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    Nextflow Multi-Omics Secondary Analysis Core
========================================================================================
*/

// Define default tracking parameters
params.reads = "data/fastq/*_{1,2}.fastq.gz"
params.outdir = "results/pipeline_output"

/*
 * PROCESS 1: Quality Control Filtering via FastQC simulation
 */
process FASTQC {
    tag "QC on ${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html_reports

    script:
    """
    echo "[Nextflow Engine] Running high-throughput quality control metrics on ${sample_id}..."
    touch "${sample_id}_1_fastqc.html"
    touch "${sample_id}_2_fastqc.html"
    """
}

/*
 * PROCESS 2: Core Transcriptome Quantification Simulation
 */
process QUANTIFICATION {
    tag "Quantifying ${sample_id}"
    publishDir "${params.outdir}/quant", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_counts.txt", emit: counts

    script:
    """
    echo "[Nextflow Engine] Aligning reads and quantifying expression maps for ${sample_id}..."
    # Simulate an expression quantification vector output
    echo "Gene_ID\t${sample_id}" > ${sample_id}_counts.txt
    echo "YAL001C\t50" >> ${sample_id}_counts.txt
    echo "YBL002W\t120" >> ${sample_id}_counts.txt
    """
}

/*
 * WORKFLOW TERMINAL CONTROL LOGIC
 */
workflow {
    log.info """
    R N A - S E Q   N E X T F L O W   P I P E L I N E
    ==================================================
    Reads Input Set : ${params.reads}
    Output Directory: ${params.outdir}
    ==================================================
    """

    // Explicitly casting strings into Nextflow file objects using file()
    channel_inputs = Channel.of(
        tuple('Sample_Control_Rep1', [file('mock_R1.fq.gz'), file('mock_R2.fq.gz')]),
        tuple('Sample_Treatment_Rep1', [file('mock_R1.fq.gz'), file('mock_R2.fq.gz')])
    )

    // 2. Execute parallel Quality Control processes
    FASTQC(channel_inputs)

    // 3. Execute Quantifications in parallel channels
    QUANTIFICATION(channel_inputs)
    
    // 4. Collect pipeline summaries
    QUANTIFICATION.out.counts.collect()
        .view { it -> "[Pipeline Complete] Matrix components successfully assembled: ${it}" }
}
