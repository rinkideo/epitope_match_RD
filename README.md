# Epitope Match

Pipeline to match reads/sequences to predefined HIV epitopes (HXB2 coordinates) using an identity threshold.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18852547.svg)](https://doi.org/10.5281/zenodo.18852547)

## Repo structure
- `epitope_match_80perc_th.py`: main script
- `slurm_epimatch_array.sh`: SLURM array runner
- `data/reference/HIV_Reference.fasta`: reference FASTA
- `data/Epitopes_HXB2_Coordinates.xlsx`: epitope coordinates table
- `data/samples.txt`: one sample ID per line (example)

## SLURM usage
Run from repo root:

sbatch --array=0-14 slurm_epimatch_array.sh

Edit paths inside `slurm_epimatch_array.sh` for:
- FASTQ directory
- sample list
- output directory

## Python usage
python epitope_match_80perc_th.py --sample_id SampleA --r1 R1.fastq.gz --r2 R2.fastq.gz --ref data/reference/HIV_Reference.fasta --epitopes data/Epitopes_HXB2_Coordinates.xlsx --output-dir results/SampleA_out
