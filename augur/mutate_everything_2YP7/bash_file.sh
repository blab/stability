#!/bin/sh
srun -n 1 -c 1 -t 72:00:00 -o output_1.txt python 1_foldx_split/run_mutator_cluster.py 1 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_2.txt python 2_foldx_split/run_mutator_cluster.py 2 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_3.txt python 3_foldx_split/run_mutator_cluster.py 3 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_4.txt python 4_foldx_split/run_mutator_cluster.py 4 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_5.txt python 5_foldx_split/run_mutator_cluster.py 5 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_6.txt python 6_foldx_split/run_mutator_cluster.py 6 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_7.txt python 7_foldx_split/run_mutator_cluster.py 7 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_8.txt python 8_foldx_split/run_mutator_cluster.py 8 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_9.txt python 9_foldx_split/run_mutator_cluster.py 9 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_10.txt python 10_foldx_split/run_mutator_cluster.py 10 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_11.txt python 11_foldx_split/run_mutator_cluster.py 11 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_12.txt python 12_foldx_split/run_mutator_cluster.py 12 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_13.txt python 13_foldx_split/run_mutator_cluster.py 13 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_14.txt python 14_foldx_split/run_mutator_cluster.py 14 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_15.txt python 15_foldx_split/run_mutator_cluster.py 15 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_16.txt python 16_foldx_split/run_mutator_cluster.py 16 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_17.txt python 17_foldx_split/run_mutator_cluster.py 17 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_18.txt python 18_foldx_split/run_mutator_cluster.py 18 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_19.txt python 19_foldx_split/run_mutator_cluster.py 19 1HA0 & 
srun -n 1 -c 1 -t 72:00:00 -o output_20.txt python 20_foldx_split/run_mutator_cluster.py 20 1HA0 & 
wait