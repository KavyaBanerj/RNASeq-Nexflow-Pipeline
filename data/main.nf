#!/usr/bin/env nextflow
 
 // All parameters for the pipeline
params.reads = "$baseDir/data/reads/*_{1,2}.fastq.gz" // Location of input reads
params.outdir = "$baseDir/results" // Output directory

params.ref_genome = "$baseDir/references/genome.fa" // Reference genome location
params.adapters = "$baseDir/config/adapters.fa" // Adapter sequences for Trimmomatic
params.star_index = "$baseDir/references/star" // STAR index directory
params.transcriptome = "$baseDir/references/transcriptome.fa" // Transcriptome file for Salmon

// Log
log.info """\
RNA-Seq Analysis Pipeline
Version: 1.0
Reads: ${params.reads}
Output directory: ${params.outdir}
"""

workflow.onComplete {
    file(params.outdir).mkdirs()
    log.info "Pipeline complete. Output files are in ${params.outdir}"
}

