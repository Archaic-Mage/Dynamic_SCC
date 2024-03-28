#!/bin/bash
#PBS -o logfile.log
#PBS -e errorfile_slash.err
#PBS -l walltime=00:05:00
#PBS -l select=2:ncpus=32

module load openmpi316

tpdir=`echo $PBS_JOBID | cut -f 1 -d .`
tempdir=$HOME/scratch/job$tpdir
mkdir -p $tempdir
cd $tempdir
cp $PBS_O_WORKDIR/build/dynamic_scc .

#Execution
mpirun -np 2 dynamic_scc < $PBS_O_WORKDIR/tests/input/i4.txt > output.txt 

mv ../job$tpdir $PBS_O_WORKDIR/.
