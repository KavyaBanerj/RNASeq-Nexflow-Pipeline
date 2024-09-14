#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// All parameters for the pipeline
params.reads =  params.reads ?: "$PWD/data/reads/*{1,2}.fastq.gz" // Location of input reads
params.outdir = params.outdir ?: "$PWD/results" // Output directory
params.genome = params.genome ?: "$PWD/data/refGenome/genome.fasta" // Path to the genome FASTA file
params.gtf = params.gtf ?: "$PWD/data/refGenome/genome.gtf" // Path to the GTF file
params.transcripts = params.transcripts ?: "$PWD/data/refGenome/transcripts.fasta"  // Path to the transcriptome FASTA file
params.skip_alignment  = params.skip_alignment ?: false
params.test_mode       = params.test_mode ?: false

// Validate parameter types and provide defaults
if (!params.containsKey('reads')) {
    error "Parameter 'reads' is required"
}

if (!file(params.genome).exists()) {
    error "Genome file not found: ${params.genome}"
}
if (!file(params.gtf).exists()) {
    error "GTF file not found: ${params.gtf}"
}
if (!file(params.transcripts).exists()) {
    error "Transcript file not found: ${params.transcripts}"
}

// Log
log.info """\
RNA-Seq Analysis Pipeline -  FastQC, TrimGalore, STARIndex, STARAlign , Salmon Quantification, and MultiQC
Reads: ${params.reads}
Output directory: ${params.outdir}
Genome: ${params.genome}
GTF: ${params.gtf}
Transcripts: ${params.transcripts}
Test mode : ${params.test_mode}
"""

// Create output directories if they do not exist
if (!file(params.outdir).exists()) {
    file(params.outdir).mkdirs()
}

// FastQC process for quality control of reads
process FastQC {
    tag "FastQC on ${sample_id}"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html"), emit: reports

    script:
    """
    fastqc -o . ${reads[0]} ${reads[1]}

     echo "FASTQC completed for ${sample_id}"
    """
}

// TrimGalore process for trimming low-quality reads
process TrimGalore {
    tag "TrimGalore on ${sample_id}"
    container 'quay.io/biocontainers/trim-galore:0.6.6--hdfd78af_1'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1_trimmed.fq.gz"), path("${sample_id}_2_trimmed.fq.gz"), emit: trimmed_reads

    script:
    """
    trim_galore --paired --quality 20 --length 36 --gzip --output_dir . ${reads[0]} ${reads[1]}
    mv ${sample_id}_1_val_1.fq.gz ${sample_id}_1_trimmed.fq.gz
    mv ${sample_id}_2_val_2.fq.gz ${sample_id}_2_trimmed.fq.gz

     echo "Trimmed Reads completed for ${sample_id}"
    """
}

// Creating Index for STAR alignment
process STARIndex {
    tag "STARIndex"
    container 'quay.io/biocontainers/star:2.7.8a--0'
    publishDir "${params.outdir}/STARindex", mode: 'copy'

    input:
    path genome_fasta
    path gtf_file

    output:
    path("STARindex"), emit: star_index

    script:
    """
    mkdir -p STARindex
    STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --genomeDir STARindex \
    --genomeFastaFiles ${genome_fasta} \
    --sjdbGTFfile ${gtf_file} \
    --sjdbOverhang 100 \
    --genomeSAindexNbases 11
    """
}

// STAR Alignment process for creating aligned bam files
process STARAlign {
    tag "STARAlign on ${sample_id}"
    container 'quay.io/biocontainers/star:2.7.8a--0'
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads), path(star_index)

    output:
    tuple val(sample_id), path("star_out/${sample_id}.Aligned.sortedByCoord.out.bam"), emit: sorted_bam
    tuple val(sample_id), path("star_out/${sample_id}.Aligned.toTranscriptome.out.bam"), emit: transcriptome_bam
    tuple val(sample_id), path("star_out/${sample_id}.Chimeric.out.junction")
    tuple val(sample_id), path("star_out/${sample_id}.Log.final.out")
    tuple val(sample_id), path("star_out/${sample_id}.Log.out")
    tuple val(sample_id), path("star_out/${sample_id}.Log.progress.out")
    tuple val(sample_id), path("star_out/${sample_id}.ReadsPerGene.out.tab")
    tuple val(sample_id), path("star_out/${sample_id}.SJ.out.tab")
    script:
    """
    mkdir -p star_out
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${star_index} \
         --readFilesIn ${trimmed_reads[0]} ${trimmed_reads[1]} \
         --readFilesCommand zcat \
         --outFileNamePrefix star_out/${sample_id}. \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts \
         --outFilterMultimapNmax 20 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.1 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --outFilterType BySJout \
         --outFilterScoreMinOverLread 0.33 \
         --outFilterMatchNminOverLread 0.33 \
         --limitSjdbInsertNsj 1200000 \
         --outSAMstrandField intronMotif \
         --outFilterIntronMotifs None \
         --alignSoftClipAtReferenceEnds Yes \
         --outSAMattrRGline ID:rg1 SM:sm1 \
         --outSAMattributes NH HI AS nM NM ch \
         --chimSegmentMin 15 \
         --chimJunctionOverhangMin 15 \
         --chimOutType Junctions WithinBAM SoftClip \
         --chimMainSegmentMultNmax 1

    echo "STAR alignment completed for ${sample_id}"
    """
}

// Salmon for quantification on reads 
process SalmonQuant {
    tag "SalmonQuant on ${sample_id}"
    container 'quay.io/biocontainers/salmon:1.10.3--hb7e2ac5_1'
    publishDir "${params.outdir}/salmon_quant", mode: 'copy'

    input:
    tuple val(sample_id), path(transcriptome_bam), path(transcript_file)

    output:
    path("${sample_id}_quant/quant.sf"), emit: quant_results
    path("${sample_id}_quant/aux_info")
    val "${sample_id}"

    script:
    """
    salmon quant -t ${transcript_file} -l A -a ${transcriptome_bam} -o ./${sample_id}_quant --gcBias

    echo "Salmon Qauntification completed for ${sample_id}"
    """
}

// MultiQC for aggregating reports
process MultiQC {
    tag "MultiQC"
    container 'quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_reports

    output:
    path("multiqc_report.html"), emit: multiqc_report

    script:
    """
    set -e
    multiqc ${qc_reports} -o ./
    """
}

// TestProcess for validation
process TestProcess {
    tag 'TestProcess'
    publishDir "${params.outdir}/test_output", mode: 'copy'

    script:
    """
    set -e
    echo "This is a test process to verify pipeline setup." > test_output.txt
    """
}

// Main workflow
workflow {
    // Define channels
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true, size: 2)
    genome_ch = Channel.value(file(params.genome))
    gtf_ch = Channel.value(file(params.gtf))
    transcript_ch = Channel.value(file(params.transcripts))

    star_index_ch = STARIndex(genome_ch, gtf_ch)

    // Call FastQC
    reads_ch.map { sample_id, files -> tuple(sample_id, files) }
            .set { fastqc_ch }

    FastQC(fastqc_ch)

    // Connect reads to Trim Galore (Note: We process reads directly without linking FastQC outputs)
    reads_ch.map { sample_id, files -> tuple(sample_id, files) }
            .set { trimmed_ch }

    // Trim Galore
    TrimGalore(trimmed_ch)

    trimmed_ch.combine(star_index_ch).map { sample_id, files, star_index -> tuple(sample_id, files, star_index) }
        .set { star_align_ch }

    // STAR Align
    star_output = STARAlign(star_align_ch)

    // Salmon Quant
    star_output.transcriptome_bam.combine(transcript_ch).map { sample_id, bam, transcript_file -> tuple(sample_id, bam, transcript_file) }
        .set { salmon_quant_ch }

    SalmonQuant(salmon_quant_ch)
}
