#!/bin/bash
#
#  This is a simple example of how to run a job under MOAB.
#
#MOAB -N hello_moab
#MOAB -j oe
#MOAB -q gunzburg_q
#
date
#
echo "qsub is running on " $PBS_O_HOST
echo "Originating queue is " $PBS_O_QUEUE
echo "Executing queue is " $PBS_QUEUE
echo "Execution mode is " $PBS_ENVIRONMENT
echo "Job identifier is " $PBS_JOBID
echo "Job name is " $PBS_JOBNAME
echo "Node file is " $PBS_NODEFILE
echo "Home directory is " $PBS_O_HOME
echo "Working directory is " $PBS_O_WORKDIR
echo "Path = " $PBS_O_PATH
#
#  Move to the directory from which you submitted the script.
#
cd $PBS_O_WORKDIR
#
gfortran hello.f90
mv a.out hello
./hello
rm hello
#
date
#
exit
