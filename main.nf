nextflow.enable.dsl=2
 process FASTQC {
    tag "Running FastQC on $sample_id"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.zip", emit: zip
    
    script:
    """
    fastqc -q ${reads[0]} ${reads[1]}
    """
}

process ALIGN {
    tag "Aligning $sample_id with HISAT2"
    
    input:
    tuple val(sample_id), path(reads)
    path genome_index
    
    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam
    
    script:
    """
    hisat2 -x genome_index -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam --no-spliced-alignment -p 2
    """
}

process SAMTOOLS {
    tag "Processing SAM/BAM for $sample_id"
    
    input:
    tuple val(sample_id), path(sam)
    
    output:
    path "${sample_id}.bam", emit: bam
    
    script:
    """
    samtools view -bS $sam | samtools sort -o ${sample_id}.bam
    """
}

process FEATURECOUNTS {
    tag "Quantifying expression matrix"
    publishDir "results", mode: 'copy'

    input:
    path bams
    path gtf

    output:
    path "sample_counts.txt"
    path "sample_counts.txt.summary"

    script:
    """
    featureCounts -a $gtf -o sample_counts.txt -p -T 2 $bams
    """
}

process MULTIQC {
    publishDir "results", mode: 'copy'

    input:
    path fastqc_files

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

process DESEQ2 {
    tag "Running differential expression"
    publishDir "results", mode: 'copy'

    input:
    path counts_file
    path r_script
 

    output:
    path "deseq2_results.csv"

    script:
    """
    Rscript ${r_script}
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads)
    genome_index_ch = Channel.fromPath("${params.genome_index}*.ht2").collect()

    FASTQC(read_pairs_ch)
    ALIGN(read_pairs_ch, genome_index_ch)
    SAMTOOLS(ALIGN.out.sam)
    FEATURECOUNTS(SAMTOOLS.out.bam.collect(), params.gtf)
    MULTIQC(FASTQC.out.zip.collect())
    
    // Pass only the primary count matrix text file into DESeq2
    DESEQ2(FEATURECOUNTS.out[0], "${baseDir}/run_deseq2.R")

}

  
 
