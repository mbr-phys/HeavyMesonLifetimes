#!/bin/bash 

#SBATCH -J FermionFlow
#SBATCH -A ACCOUNT
#SBATCH --time=4:00:00
#SBATCH --partition=gpu
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err
#SBATCH --no-requeue

Lattice="32.32.32.64"
MPIGrid="1.1.2.4"

machine=MACHINE
input=input_fermion

echo "--------------------------------"
echo "SLURM job running on: `hostname`"
echo "in directory:         `pwd`"
echo "SLURM jobid:          ${SLURM_JOB_ID}"
echo "SLURM #nodes:         ${SLURM_NNODES}"
echo "SLURM tasks per node: ${SLURM_TASKS_PER_NODE}"
echo "SLuRM #ntasks:        ${SLURM_NTASKS}" 
echo "Nodefile: ${SLURM_JOB_NODELIST}"
echo "--------------------------------"

echo "# time-start "`date`
TotalTic=`date +%s`

##################################################
##### SET ENV ARGS AND LOAD MODULES ETC HERE #####
##################################################

HADRONS=/path/to/HeavyMesonLifetimes/build/FermionFlow_Lifetimes
PARAMS=" $input.xml --accelerator-threads 8 --grid ${Lattice} --mpi ${MPIGrid} --comms-concurrent --comms-overlap --shm 2048 --shm-mpi 1 --dslash-unroll --decomposition"

echo "srun -n ${SLURM_NTASKS} --cpu-bind=${CPU_BIND} ${SLURM_SUBMIT_DIR}/select_gpu ${HADRONS} ${PARAMS}"
mpirun -np ${SLURM_NTASKS} -x LD_LIBRARY_PATH --bind-to none \
    ./gpu-mpi-wrapper.sh \
    ${HADRONS} ${PARAMS} 

############################################################################################# 
############################################################################################# 

TotalToc=`date +%s`
echo "# time-finish "`date`

TotalTime=$(( $TotalToc - $TotalTic ))
TotalHours=`echo "$TotalTime / 3600" | bc -l`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Total time  $TotalTime [sec] = $TotalHours [h]"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
