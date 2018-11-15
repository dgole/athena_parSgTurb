#!/bin/bash

#SBATCH -J midTest
##SBATCH --time=48:00:00
#SBATCH --time=1:00:00
#SBATCH -N 1
#SBATCH --ntasks 64
#SBATCH -o out_%j.txt
#SBATCH -e error_%j.txt
##SBATCH -p normal
#SBATCH -p development

# The following commands will be executed when this script is run.
module load intel impi fftw3 
ibrun -np 64 ./athena -i athinput.par_strat3d_turb 


#
