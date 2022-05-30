#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=sphere_sim
# Archivo de salida
#SBATCH --output=sphere.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=pzurita@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 20 python3 sphere_main.py

