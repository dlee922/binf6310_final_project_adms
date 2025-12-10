#!/bin/bash
#SBATCH --job-name=tem1_expression
#SBATCH --partition=courses
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=expression_model_%j.out
#SBATCH --error=expression_model_%j.err

module load R

cd /courses/BINF6310.202610/students/lo.mich/final_group_project

Rscript run_expression_model.R
