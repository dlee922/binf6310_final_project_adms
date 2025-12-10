#!/bin/bash
#SBATCH --job-name=save_mic_model
#SBATCH --partition=courses
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=mic_safe_%j.out
#SBATCH --error=mic_safe_%j.err

module load R

cd /courses/BINF6310.202610/students/lo.mich/final_group_project

Rscript modelMIC_ROBUST.R

