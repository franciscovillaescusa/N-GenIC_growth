#!/bin/bash
#SBATCH -p general
#SBATCH -J Run_0.0eV 
#SBATCH -n 264
#SBATCH --exclusive
#SBATCH -o OUTPUT.o%j
#SBATCH -e OUTPUT.e%j
#SBATCH --mail-user=villaescusa.francisco@gmail.com
#SBATCH --mail-type=ALL

mpirun ../n-genic_growth/N-GenIC N-GenIC.param >> logIC
mpirun ../g3_yb_hz/P-Gadget3 ics.param >> logfile
