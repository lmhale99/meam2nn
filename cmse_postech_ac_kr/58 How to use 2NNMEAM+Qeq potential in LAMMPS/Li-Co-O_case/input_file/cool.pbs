#!/bin/sh

### Job name ##########
#PBS -N lammps_Ex


### Queue name #########
##PBS -q  default

### Number of nodes ####
#PBS -l nodes=1:ppn=8

### This job's working directory ##########
cd $PBS_O_WORKDIR    
NPROCS=$(wc -l < $PBS_NODEFILE)

###########################################
# Run your executable
###########################################

mpirun -np $NPROCS lmp_mkl -i in.cool -l log.cool >& out

