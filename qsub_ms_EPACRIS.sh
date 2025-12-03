#!/bin/bash
#PBS -N EC-test
#PBS -q array-sn
##PBS -I -q interactive-jpl
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=24:00:00
#PBS -k oed
##PBS -o job_stdout.txt
##PBS -e job_stderr.txt
#PBS -W group_list=renyu-group
#PBS -v summary=true
 
### Load modules into your environment
# module load intel/compiler/64
# module load intel/mkl/64
 
### Run executable
echo "$PBS_JOBID"
#cd /scratch_lg/renyu-group/mscheucher/EPACRIS-1/
cd /scratch_lg_edge/renyu-group/mscheucher/Epacris/
./Epacris $PBS_JOBID $PBS_JOBNAME
