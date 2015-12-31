#!/bin/sh
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --output=batch_output.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=COMPLETED
#SBATCH --mail-user=chacalle@uw.edu
#SBATCH --mail-user=chacalle@uw.edu

srun python mutator_stability.py