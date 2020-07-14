#!/bin/bash
#SBATCH --mem=50G
#SBATCH --time=2:30:0
#SBATCH --error=logs/err/get_gv_pairs.err
#SBATCH --output=logs/out/get_gv_pairs.out

module load R
Rscript ../code/get_gv_pairs.R
