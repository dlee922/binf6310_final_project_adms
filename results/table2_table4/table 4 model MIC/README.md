These three files form the complete workflow for the MIC (minimum inhibitory concentration) 
model that generates Table 4:


## Scripts
- modelMIC_ROBUST.R - The core model script that Prepares data (292 isolates final dataset), Fits MCMCglmm ordinal model, Generates Table 4 results
- save_output_as_it_goes_submit_mic.slurm - SLURM submission (24hr, 32GB memory)

## Output file
- Output folder (graphs)
- mic_safe_3042794.out - table 4
