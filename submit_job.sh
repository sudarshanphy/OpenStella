#!/bin/bash

#SBATCH -A ast137
#SBATCH -q debug
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -J 1_openstella_new3_s15_1718 
#SBATCH -o 1_openstella_new3_s15_1718.%j.out
#SBATCH -e 1_openstella_new3_s15_1718.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sneopane@vols.utk.edu 

module load PrgEnv-gnu 
date
echo "OpenStella run for s15 1718 flashx frame"
echo "cells with velx > 1.0e6 selected"
echo "inner 56 cells are mapped as it is"
echo "remianing celles are combined from 2 to 1"
echo "--------------------------------"
HOMEStella=$PWD
export HOMEStella

cd $HOMEStella/run/strad || exit 1

#rm -f st.log

echo "Starting STELLA in:"
pwd
date

#srun -N 1 -n 1 ./xstella6y12m.exe > st.log 2>&1
srun -N 1 -n 1 ./xstella6y12m.exe

echo "STELLA finished"

date

