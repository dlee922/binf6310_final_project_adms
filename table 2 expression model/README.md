This directory contains the computational workflow for reproducing 
Table 2 from Matlock et al. (Nature Communications, 2025), which 
presents a Bayesian phylogenetic mixed-effects model of TEM-1 β-lactamase 
expression levels in E. coli as a function of promoter variants and gene copy 
number.




## Files Overview
- submit_expression.slurm - SLURM job submission script for HPC execution
- run_expression_model.R - Wrapper script that handles package installation and runs the analysis
- modelExpression.R - Main R script containing all data processing and model fitting code
- expression_model_3016329.out - HPC output log file showing job execution and convergence diagnostics

## need data from ./data/
- supplementary_data_1.tsv
- supplementary_data_2.tsv
- ml_tree.txt
- supplementary_data_3.xlsx
