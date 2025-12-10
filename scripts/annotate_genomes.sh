#!/bin/bash
#SBATCH --job-name=annotate_ecoli
#SBATCH --output=logs/annotate_%A_%a.out
#SBATCH --error=logs/annotate_%A_%a.err
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-50

echo "****************************************************"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Run on host: $(hostname)"
echo "Started at: $(date)"
echo "****************************************************"

# Load conda and activate environment
module load miniconda3
source activate ecoli_final

# Project directory
PROJECT_DIR=/scratch/evangelatos.s/final_project
cd $PROJECT_DIR

# Get the genome file for this array task
GENOME_LIST=(assemblies/*.fna)
GENOME=${GENOME_LIST[$SLURM_ARRAY_TASK_ID - 1]}
SAMPLE=$(basename $GENOME .fna)

echo "Processing sample: $SAMPLE"
echo "Genome file: $GENOME"

# Create output directories
mkdir -p annotations/${SAMPLE}
mkdir -p amrfinder_output
mkdir -p mlst_output

echo "****************************************************"
echo "Running Prokka annotation at: $(date)"

# Run Prokka with parameters from paper
prokka \
    --outdir annotations/${SAMPLE} \
    --prefix ${SAMPLE} \
    --cpus 8 \
    --centre X \
    --compliant \
    --genus Escherichia \
    --species coli \
    --force \
    $GENOME

echo "Prokka complete at: $(date)"

echo "****************************************************"
echo "Running AMRFinder at: $(date)"

# Run AMRFinder to find resistance genes
amrfinder \
    -n $GENOME \
    -o amrfinder_output/${SAMPLE}.tsv \
    --plus \
    --organism Escherichia \
    --threads 8

echo "AMRFinder complete at: $(date)"

echo "****************************************************"
echo "Running MLST at: $(date)"

# Run MLST to verify sequence type
mlst $GENOME > mlst_output/${SAMPLE}.tsv

echo "MLST complete at: $(date)"

echo "****************************************************"
echo "Finished sample $SAMPLE at: $(date)"
echo "Output files:"
echo "  - Prokka: annotations/${SAMPLE}/"
echo "  - AMRFinder: amrfinder_output/${SAMPLE}.tsv"
echo "  - MLST: mlst_output/${SAMPLE}.tsv"
echo "****************************************************"
