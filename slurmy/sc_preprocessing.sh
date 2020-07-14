#!/bin/bash
#SBATCH --mem=50G
#SBATCH --time=3:0:0
#SBATCH --error=logs/err/sc_pre.err
#SBATCH --output=logs/out/sc_pre.out

module load R
Rscript ../code/sc_preprocessing.R
