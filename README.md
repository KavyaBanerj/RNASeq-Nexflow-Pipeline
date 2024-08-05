
# RNA-Seq QC and Quantification Pipeline

This repository contains a Nextflow pipeline for RNA-Seq QC and quantification on pair-ended reads. The pipeline includes quality control, read trimming, alignment, and quantification steps.

## Pipeline Overview

The pipeline performs the following steps:
1. **Quality Control** using FastQC
2. **Read Trimming** using TrimGalore
3. **Genome Indexing** using STAR
4. **Read Alignment** using STAR
5. **Quantification** using Salmon

## Prerequisites
- **Nextflow**: Ensure you have Nextflow installed. Follow the installation instructions [here](https://www.nextflow.io/docs/latest/getstarted.html).
- **Docker:** The pipeline uses containers for reproducibility. Install Docker [here](https://docs.docker.com/get-docker/)
- Recommended to load the pipeline on a Gitpod work environment.

## Input Files

The pipeline requires the following input files:
- Paired-end RNA-Seq reads in FASTQ format
- Reference genome in FASTA format
- Gene annotation file in GTF format
- Transcriptome FASTA file

## Usage

### Parameters

- `reads`: Location of input reads (default: `$PWD/data/reads/*{1,2}.fastq.gz`)
- `outdir`: Output directory (default: `$PWD/results`)
- `genome`: Path to the genome FASTA file (default: `$PWD/data/refGenome/genome.fasta`)
- `gtf`: Path to the GTF file (default: `$PWD/data/refGenome/genome.gtf`)
- `transcripts`: Path to the transcriptome FASTA file (default: `$PWD/data/refGenome/transcript.fasta`)

### Running the Pipeline

To run the pipeline with default parameters, use the following command:

```bash
nextflow run main.nf
```

This will use the default paths specified in the pipeline script for the input files and output directory.

### Example Command with Custom Parameters

To specify custom input files and output directory, use the following command:

```bash
nextflow run main.nf --reads <path_to_reads> --genome <path_to_genome> --gtf <path_to_gtf> --transcripts <path_to_transcripts> --outdir <output_directory>
```

## Output

The pipeline generates the following outputs:
- **FastQC Reports**: Quality control reports for raw reads (`<outdir>/qc/`)
- **Trimmed Reads**: Trimmed read files (`<outdir>/trimmed/`)
- **STAR Index**: STAR index files (`<outdir>/STARindex/`)
- **Aligned BAM files**: Aligned and sorted BAM files (`<outdir>/aligned/`)
- **Salmon Quantification**: Quantification results (`<outdir>/salmon_quant/`)

## Pipeline Processes

### FastQC

Performs quality control on raw reads.

### TrimGalore

Trims low-quality bases and adapters from reads.

### STARIndex

Creates a STAR genome index from the reference genome and GTF file.

### STARAlign

Aligns the trimmed reads to the reference genome using STAR.

### SalmonQuant

Quantifies transcript abundance using Salmon using a transcript FASTA file.

## Dependencies

The pipeline uses the following software containers from [Biocontainers] (https://biocontainers.pro/) :

- **FastQC**: `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0`
- **TrimGalore**: `quay.io/biocontainers/trim-galore:0.6.6--hdfd78af_1`
- **STAR**: `quay.io/biocontainers/star:2.7.8a--0`
- **Salmon**: `quay.io/biocontainers/salmon:1.10.3--hb7e2ac5_1`

Ensure these containers are available in your Docker environment.