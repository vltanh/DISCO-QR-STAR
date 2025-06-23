#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --output=slurm_output/slurm-%j.out
#SBATCH --job-name="agg_result"
#SBATCH --partition=secondary
#SBATCH --mem=32G

/usr/bin/time -v python agg_result.py --input output/trees/ --output ncd.csv --reference data/trees/