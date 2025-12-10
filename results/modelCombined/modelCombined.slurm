#!/bin/bash
#SBATCH --job-name=model_combined      # Name of the job
#SBATCH --partition=courses                # Who to bill for the job
#SBATCH -N 1                              # How many nodes do you need. When in doubt, use 1
#SBATCH -c 8                              # How many "threads" do you need for your job
#SBATCH --mem 48G                         # How much memory
#SBATCH -t 24:00:00                        # How long
#SBATCH --out=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err

# Load env
module load miniconda3
conda activate tem_env
module load R/4.4.1

# Point R to your library
export R_LIBS_USER="/home/desai.ara/R/x86_64-pc-linux-gnu-library/4.4"

# Run your R script
Rscript /home/desai.ara/binf6310_final_project_adms/scripts/modelCombined.R