# slurm_epimatch_array.sh
#
# Generic SLURM array runner for the epitope-matching pipeline.
# Notes:
# - This script expects paired FASTQs named:
#     <SAMPLE_ID>_R1_reads.fastq.gz and <SAMPLE_ID>_R2_reads.fastq.gz

#!/bin/bash
#SBATCH --job-name=epimatch
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=5-00:00:00
#SBATCH --array=0-14
#SBATCH --output= %x_%A_%a.out
#SBATCH --error= %x_%A_%a.err

set -euo pipefail

# Run from repo root (recommended):
WORKDIR="$(pwd)"
cd $WORKDIR

#cluster-specific environment setup
module load gcc/14.2.0
module load python/3.13.1
source "${HOME}/biopython_env/bin/activate"

SAMPLE_LIST= data/samples.txt
SAMPLE_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SAMPLE_LIST)

if [ -z "$SAMPLE_ID" ]; then
  echo "No sample found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID. Exiting."
  exit 0
fi


FASTQ_DIR=data/fastq_reads
R1=${FASTQ_DIR}/${SAMPLE_ID}_R1_reads.fastq.gz
R2=${FASTQ_DIR}/${SAMPLE_ID}_R2_reads.fastq.gz
OUTDIR= results/${SAMPLE_ID}_out
HIV_REF=data/reference/HIV_Reference.fasta
EPITOPE_FILE=data/Epitopes_HXB2_Coordinates.xlsx


##Remove old output if it exists
if [ -d "$OUTDIR" ]; then
 echo "Removing existing output directory: $OUTDIR"
 rm -rf "$OUTDIR"
fi

mkdir -p $OUTDIR


python epitope_match_80perc_th.py --sample_id "$SAMPLE_ID" --r1 "$R1" --r2 "$R2" --ref $HIV_REF --epitopes $EPITOPE_FILE --output-dir "$OUTDIR"



