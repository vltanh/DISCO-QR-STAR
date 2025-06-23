#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --output=slurm_output/backup/slurm-%j.out
#SBATCH --job-name="project"
#SBATCH --partition=eng-instruction
#SBATCH --account=25sp-cs581a-eng
#SBATCH --mem=8G

zip -r discoqr_discoqrstar.zip output/trees/
# unzip discoqr_discoqrstar.zip