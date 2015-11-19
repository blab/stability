#!/bin/sh
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --output=terminal_output.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=COMPLETED
#SBATCH --mail-user=charltonsee@ymail.com

python foldx_cluster_runpipeline.py 4 0