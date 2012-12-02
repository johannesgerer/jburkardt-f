#!/bin/bash
#
#PBS -lwalltime=00:00:30
#PBS -lnodes=2:ppn=2
#PBS -W group_list=tcf_user
#PBS -q production_q
#PBS -A admn0000

NUM_NODES=`/bin/cat $PBS_NODEFILE | /usr/bin/wc -l | /usr/bin/sed "s/ //g"`

cd $PBS_O_WORKDIR

export PATH=/nfs/software/bin:$PATH

echo "Run program."

jmdrun -np $NUM_NODES \
  -hostfile $PBS_NODEFILE \
  ./random_mpi &> random_sysx_output.txt

echo "Clean up."

date

exit;

