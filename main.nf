#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// All parameters for the pipeline
params.reads = "$PWD/data/reads/ENCSR000COQ1_{1,2}.fastq.gz" // Location of input reads
params.outdir = "$PWD/results" // Output directory
params.genome = "$PWD/data/refGenome/genome.fasta" // Path to the genome FASTA file
params.gtf = "$PWD/data/refGenome/genome.gtf" // Path to the GTF file
params.transcripts = "$PWD/data/refGenome/transcript.fasta" // Path to the transcriptome FASTA file

// Log
log.info """\
RNA-Seq Analysis Pipeline - STARIndex, STARAlign, and Salmon
Reads: ${params.reads}
Output directory: ${params.outdir}
Genome: ${params.genome}
GTF: ${params.gtf}
Transcripts: ${params.transcripts}
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
    path("*_fastqc.zip"), emit: reports

    script:
    """
    fastqc -o . ${reads[0]} ${reads[1]}
    """
}

process TrimGalore {
    tag "TrimGalore on ${sample_id}"
    container 'quay.io/biocontainers/trim-galore:0.6.6--hdfd78af_1'
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1_trimmed.fq.gz"), path("${sample_id}_2_trimmed.fq.gz")

    script:
    """
    trim_galore --paired --quality 20 --length 36 --gzip --output_dir . ${reads[0]} ${reads[1]}
    mv ${sample_id}_1_val_1.fq.gz ${sample_id}_1_trimmed.fq.gz
    mv ${sample_id}_2_val_2.fq.gz ${sample_id}_2_trimmed.fq.gz
    """
}

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
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir STARindex --genomeFastaFiles ${genome_fasta} --sjdbGTFfile ${gtf_file} --sjdbOverhang 100 --genomeSAindexNbases 11
    """
}

process STARAlign {
    tag "STARAlign on ${sample_id}"
    memory '64 GB'
    container 'quay.io/biocontainers/star:2.7.8a--0'
    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads), path(star_index)

    output:
    tuple val(sample_id), path("star_out/${sample_id}.Aligned.sortedByCoord.out.bam"), emit: sorted_bam
    tuple val(sample_id), path("star_out/${sample_id}.Aligned.toTranscriptome.out.bam"), emit: transcriptome_bam

    script:
    """
    mkdir -p star_out
    STAR --runThreadN 16 \
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
    ls -lh star_out
    """
}

process SalmonQuant {
    tag "SalmonQuant on ${sample_id}"
    container 'quay.io/biocontainers/salmon:1.5.1--h84f40af_0'
    publishDir "${params.outdir}/salmon_quant", mode: 'copy'

    input:
    tuple val(sample_id), path(transcriptome_bam), val(transcript_file)

    output:
    path("${sample_id}_quant/quant.sf")
    path("${sample_id}_quant/aux_info")
    val "${sample_id}"

    script:
    """
    echo "Transcript file path: ${transcript_file}"
    echo "Checking if transcript file exists..."
    abs_path=\$(readlink -f ${transcript_file})
    ls -l \${abs_path}

    if [ ! -f \${abs_path} ]; then
        echo "Transcript file does not exist at \${abs_path}"
        exit 1
    fi

    salmon quant -t \${abs_path} -l A -a ${transcriptome_bam} -o ./${sample_id}_quant --gcBias
    """
}

// Define the workflow
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

    // Connect FastQC to Trim Galore (Note: We process reads directly without linking FastQC outputs)
    reads_ch.map { sample_id, files -> tuple(sample_id, files) }
            .set { trimmed_ch }

    TrimGalore(trimmed_ch)

    trimmed_ch.combine(star_index_ch).map { sample_id, files, star_index -> tuple(sample_id, files, star_index) }
        .set { star_align_ch }

    star_output = STARAlign(star_align_ch)

    star_output.transcriptome_bam.combine(transcript_ch).map { sample_id, bam, transcript_file -> tuple(sample_id, bam, transcript_file) }
        .set { salmon_quant_ch }

    SalmonQuant(salmon_quant_ch)
}