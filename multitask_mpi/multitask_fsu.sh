#!/bin/bash
#
#  Name this job.
#
#MOAB -N multitask
#
#  Run this job in the "classroom" queue.
#
#MOAB -q classroom
#
#  Request 3 processors.
#
#MOAB -l nodes=1:ppn=3
#
#  Maximum wallclock time :( Hours : Minutes : Seconds )
#
#MOAB -l walltime=00:02:00
#
#  Join the two output files (standard output and standard error) into one.
#  Based on the job name, this means the single output file will be called
#  "multitask.oJOBID" where JOBID is a job id number assigned when the job is
#  submitted.
#
#MOAB -j oe
#
#  This command is required to set up the Gnu version of OpenMPI.
#
module load gnu-openmpi
#
#  This command moves from your home directory to the directory
#  from which this script was submitted.
#
cd $PBS_O_WORKDIR 
#
#  Compile and load.
#  Wise people do this BEFORE the run!
#
mpif90 multitask_mpi.f90
#
if [ $? -ne 0 ]; then
  echo "Errors compiling multitask_mpi.f90"
  exit
fi
#
#  Rename the executable.
#
mv a.out multitask
#
#  Ask MPI to use 3 processes to run your program.
#
mpirun -np 3 ./multitask > multitask_fsu_output.txt
#
#  Clean up.
#
rm multitask
#
echo "Program output written to multitask_fsu_output.txt"

