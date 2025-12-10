#!/bin/bash
#SBATCH --job-name=iqtree_auto
#SBATCH --output=logs/iqtree_auto_%j.out
#SBATCH --error=logs/iqtree_auto_%j.err
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00

module load miniconda3
source activate ecoli_final

WORK_DIR=/scratch/evangelatos.s/final_project/core_genes_v2
ALN_FILE=${WORK_DIR}/concatenated_alignment.fasta
OUT_DIR=${WORK_DIR}/iqtree_output

# Clean and recreate output directory
rm -rf $OUT_DIR
mkdir -p $OUT_DIR
cd $OUT_DIR

echo "Building tree at: $(date)"
echo "Alignment: $ALN_FILE"

# Try simpler approach - let IQ-TREE auto-select model
iqtree -s $ALN_FILE \
       -m MFP \
       -B 1000 \
       --prefix ecoli_50genomes \
       -T 16

echo "Complete at: $(date)"
ls -lh
