# Automated NGS Quality Control Pipeline

A production-ready, highly reproducible bioinformatics pipeline built with **Nextflow (DSL2)** to automate quality control checks on Next-Generation Sequencing (NGS) data.

## 🚀 Features
* **Parallel Processing:** Efficiently scales to handle paired-end sequencing reads simultaneously.
* **Automated FastQC Analysis:** Generates individual sequence quality metrics and comprehensive HTML reports.
* **Integrated MultiQC Dashboard:** Aggregates separate analysis logs into a single, interactive, summary report.
* **Resumable Core Logic:** Leverages Nextflow's native caching mechanism to skip previously computed steps, drastically reducing compute overhead.
* **Splice-Aware Coordinate Alignment:** Integrates `HISAT2` to map high-quality reads directly against an indexed reference genome.
* **Efficient Output Compression:** Leverages `Samtools` to stream, sort, and index bulky alignment data into highly compressed binary `BAM` formats.
---

## 🛠️ Pipeline Architecture
The pipeline uses Nextflow channels to automatically stream paired-end data through decoupled processing containers:



1. **Inputs:** Raw paired-end `fastq.gz` files.
2. **Process 1 (FASTQC):** Evaluates sequence quality scores, adapter contamination, and GC distribution across individual files.
3. **Process 2 (MULTIQC):** Compiles standalone metrics into a single unified client deliverable (`multiqc_report.html`).
4. **Process 3 (ALIGN):** Maps paired-end reads to a reference genome, producing dynamic text alignment diagnostics (`.summary.txt`).
5. **Process 4 (SAMTOOLS):** Converts text-heavy SAM streams into a sorted, indexed binary configuration (`.sorted.bam` and `.sorted.bam.bai`) for downstream processing.
---

## 💻 How to Run Locally

### Prerequisites
Ensure your environment (WSL2, Linux, or HPC cluster) has the following dependencies installed:
* **Java** (v11 or later)
* **Nextflow**
* **FastQC**
* **MultiQC**
* **HISAT2**
* **Samtools** 
### Execution
Clone this repository, navigate to the folder, and run the pipeline using:
```bash
nextflow run main.nf

nextflow run main.nf -resume
