#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: bash run_rnaseq.sh sample_info.csv raw_dir genome.fa genes.gtf outdir"
  exit 1
fi

META="$1"
RAW_DIR="$2"
GENOME="$3"
GTF="$4"
OUTDIR="$5"

THREADS="${THREADS:-8}"
mkdir -p "$OUTDIR"/{qc/fastp,qc/fastqc,qc/multiqc,trimmed,index,bam,counts,logs}

echo "[1/8] Validate inputs"
[[ -f "$META" ]] || { echo "Missing metadata: $META"; exit 1; }
[[ -d "$RAW_DIR" ]] || { echo "Missing raw data dir: $RAW_DIR"; exit 1; }
[[ -f "$GENOME" ]] || { echo "Missing genome fasta: $GENOME"; exit 1; }
[[ -f "$GTF" ]] || { echo "Missing GTF: $GTF"; exit 1; }

python scripts/python/validate_sample_info.py "$META"

echo "[2/8] Run fastp"
tail -n +2 "$META" | while IFS=, read -r sample group fastq1 fastq2 batch; do
  echo "Processing $sample"
  fastp \
    -i "$fastq1" \
    -I "$fastq2" \
    -o "$OUTDIR/trimmed/${sample}_R1.clean.fastq.gz" \
    -O "$OUTDIR/trimmed/${sample}_R2.clean.fastq.gz" \
    -h "$OUTDIR/qc/fastp/${sample}.fastp.html" \
    -j "$OUTDIR/qc/fastp/${sample}.fastp.json" \
    -w "$THREADS" \
    2> "$OUTDIR/logs/${sample}.fastp.log"
done

echo "[3/8] Run FastQC"
fastqc -t "$THREADS" -o "$OUTDIR/qc/fastqc" "$OUTDIR"/trimmed/*.clean.fastq.gz || true

echo "[4/8] Build HISAT2 index if needed"
INDEX_PREFIX="$OUTDIR/index/genome"
if [[ ! -f "${INDEX_PREFIX}.1.ht2" && ! -f "${INDEX_PREFIX}.1.ht2l" ]]; then
  hisat2-build -p "$THREADS" "$GENOME" "$INDEX_PREFIX"
fi

echo "[5/8] Align reads"
tail -n +2 "$META" | while IFS=, read -r sample group fastq1 fastq2 batch; do
  hisat2 -p "$THREADS" -x "$INDEX_PREFIX" \
    -1 "$OUTDIR/trimmed/${sample}_R1.clean.fastq.gz" \
    -2 "$OUTDIR/trimmed/${sample}_R2.clean.fastq.gz" \
    2> "$OUTDIR/logs/${sample}.hisat2.log" \
  | samtools view -@ "$THREADS" -bS - \
  | samtools sort -@ "$THREADS" -o "$OUTDIR/bam/${sample}.sorted.bam"

  samtools index "$OUTDIR/bam/${sample}.sorted.bam"
done

echo "[6/8] featureCounts"
featureCounts -T "$THREADS" -p -s 0 -a "$GTF" \
  -o "$OUTDIR/counts/gene_counts.txt" \
  "$OUTDIR"/bam/*.sorted.bam \
  2> "$OUTDIR/logs/featureCounts.log"

echo "[7/8] MultiQC"
multiqc "$OUTDIR" -o "$OUTDIR/qc/multiqc" || true

echo "[8/8] DESeq2 analysis"
Rscript .agents/skills/rnaseq/scripts/deseq2_basic.R \
  "$OUTDIR/counts/gene_counts.txt" \
  "$META" \
  "$OUTDIR/deseq2"

echo "RNA-seq workflow completed: $OUTDIR"
