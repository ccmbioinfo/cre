#!/bin/bash
#SBATCH -A slurm
#SBATCH --partition=all
#SBATCH --job-name=crg2
#SBATCH --time=5:00:00
#SBATCH --ntasks-per-node=4
#SBATCH --mem=20G
#SBATCH --output=%x-%j.out

module load python

python3 generate_db_from_reports.py
