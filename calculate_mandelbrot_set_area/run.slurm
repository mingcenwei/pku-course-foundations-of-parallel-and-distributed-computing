#!/bin/bash

#SBATCH --output='job_%j_%N_results.txt'
#SBATCH --partition='cpu'
#SBATCH --job-name='Mingcen-Wei_Homework-3'
#SBATCH --cpus-per-task=16
#SBATCH --time='00:01:00'

export OMP_NUM_THREADS=8
./calculate_mandelbrot_set_area
