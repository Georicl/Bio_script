# BioScript: Integrated Bioinformatics Toolbox

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Bioinformatics](https://img.shields.io/badge/field-Bioinformatics-orange.svg)](https://en.wikipedia.org/wiki/Bioinformatics)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**BioScript** (also known as `bio_tools`) is a comprehensive toolkit designed for bioinformatics analysis, focusing on sequence manipulation, RNA-seq pipelines, and sex-linked marker identification (CQ analysis). It integrates industry-standard tools into a streamlined CLI interface.

---

## 🌟 Key Features

### 1. Sequence Tools (`seqtools`)
- **Quality Assessment**: Calculate N50/L50 for genome assemblies (supports Trinity format).
- **Manipulation**: Sort sequences by length, shuffle FASTA files, and random sequence extraction.
- **Transcriptomics**: Extract the longest transcript/CDS from GFF3 and FASTA files.
- **Filtering**: Identify and filter BLAST results based on identity thresholds.

### 2. RNA-Seq Pipeline
- Full upstream to downstream workflow:
  1. HISAT2 indexing and alignment.
  2. `featureCounts` for quantification.
  3. Differential Expression (DE) analysis using **DESeq2** (with biological replicates) or **edgeR** (without replicates).
- Automated matrix generation and TMM normalization.

### 3. Sex-Linked Marker Analysis (CQ Tools)
- **`cqmapping`**: High-throughput BWA-based alignment optimized for Sex-linked Coverage Quotient (CQ) analysis.
- **`cqtools`**: Calculate CQ values by comparing female and male coverage to identify sex-specific genomic regions or scaffolds.

---

## 🛠 Prerequisites & Installation

### Software Dependencies
- **Bio-alignment**: `bwa`, `hisat2`, `samtools`, `bedtools`, `blast+`
- **R Environment**: `limma`, `edgeR`, `DESeq2`, `featureCounts` (subread package)

### Setup with `uv` (Recommended)
This project uses [uv](https://docs.astral.sh/uv/) for modern Python package management.

```bash
# 1. Clone the repository
git clone https://github.com/YourUsername/bioscript.git
cd bioscript

# 2. Sync dependencies and create virtual environment automatically
uv sync

# 3. Use the CLI tool directly
uv run bio-tools --help
```

### Traditional Setup
```bash
pip install .
bio-tools --help
```

---

## 📖 Usage Guide

BioScript uses a sub-command structure: `python main.py <command> [options]`

### Sequence Manipulation
```bash
# Calculate N50 and GC content
python main.py n50 input.fasta output_dir/

# Extract longest CDS per gene using GFF3
python main.py longest cds.fasta output_dir/ annotation.gff3

# Sort sequences (shortest to longest, use -r for reverse)
python main.py sort input.fasta output_dir/ -r
```

### RNA-Seq Workflow
The pipeline automates indexing, mapping, and DE analysis.
```bash
python main.py RNAseq \
    --gene_fasta ref.fa \
    --gtf_path ref.gtf \
    --sample_path samples.txt \
    --cpu 16 \
    --method DESeq2 \
    --data_path ./raw_data \
    --out_path ./results
```

### CQ Analysis Workflow
Used for identifying sex-linked scaffolds.
```bash
# Step 1: Alignment for both sexes (Female and Male)
python main.py cqmapping --fasta ref.fa --pair1 F_1.fq --pair2 F_2.fq -o Female/
python main.py cqmapping --fasta ref.fa --pair1 M_1.fq --pair2 M_2.fq -o Male/

# Step 2: Calculate CQ Values
python main.py cqtools \
    --f_bam Female/output.sort.bam \
    --m_bam Male/output.sort.bam \
    --fasta ref.fa \
    --output ./cq_results \
    --cq_value 0.3
```

---

## 📂 Project Structure

```text
bioscript/
├── main.py                # Command-line entry point
├── CQ_mapping/            # BWA alignment wrappers
├── CQ_tools/              # Coverage and CQ calculation logic
├── RNA_seq/               # Shell & Support scripts (R, Perl)
│   ├── RNA_seq.sh         # Main pipeline script
│   └── support_script/    # Quantification & DE analysis
└── seqtools/              # FASTA/GFF/BLAST processing modules
```

---

## ✉️ Contact & Support

- **Author**: Xiang Yang (向旸)
- **Email**: [Georicl@outlook.com](mailto:Georicl@outlook.com)
- **Issues**: Please report bugs via the [GitHub Issues](https://github.com/YourUsername/bioscript/issues) page.

---

### Citation
If you use BioScript in your research, please cite this repository:
> Xiang, Y. (2026). BioScript: An integrated pipeline for pangenome and transcriptomic analysis. GitHub: https://github.com/YourUsername/bioscript
