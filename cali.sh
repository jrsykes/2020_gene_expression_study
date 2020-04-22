#!/bin/bash
#SBATCH --partition=long
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --output=/home/sykesj/scripts/StdOut/R-%x.%j.out
#SBATCH --error=/home/sykesj/scripts/StdOut/R-%x.%j.err

python3 /home/sykesj/scripts/2020_gene_expression_study/cali.py /scratch/projects/sykesj/CaliWD/in/* /home/sykesj/dat/SRA_list_refined.csv /scratch/projects/sykesj/CaliWD/out/