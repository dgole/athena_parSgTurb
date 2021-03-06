#!/bin/bash

#SBATCH -J fp10
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks 64
#SBATCH -o out_%j.txt
#SBATCH -e error_%j.txt
#SBATCH -p normal

# The following commands will be executed when this script is run.
module load intel impi fftw3
ibrun -np 64 ./athena -i athinput.par_strat3d_turb


#
