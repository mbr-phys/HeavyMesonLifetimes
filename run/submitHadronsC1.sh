#!/bin/bash 
#SBATCH -p debug
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --time=00:59:50
#SBATCH -x hpc-node[127,179]

Lattice="24.24.24.64"
MPIGrid="1.2.2.4"
RUN=gcc

machine=omni
input=input_fermion

exec > ${SLURM_SUBMIT_DIR}/${SLURM_JOB_NAME}-${SLURM_JOB_ID}.${machine}.out 2>&1

echo "--------------------------------"
echo "SLURM job running on: `hostname`"
echo "in directory:         `pwd`"
echo "SLURM jobid:          ${SLURM_JOB_ID}"
echo "SLURM #nodes:         ${SLURM_NNODES}"
echo "SLURM tasks per node: ${SLURM_TASKS_PER_NODE}"
echo "SLURM #ntasks:        ${SLURM_NTASKS}" 
echo "Nodefile: ${SLURM_JOB_NODELIST}"
echo "--------------------------------"
lscpu

echo "# time-start "`date`
TotalTic=`date +%s`

module load gnu12/12.3.0
echo "Running Grid AVX2 with MPI on omni (gcc12.3.0/openmpi4.1.6)"
MPIRUN="mpirun "
HADRONS=<software_path>/TwoPtFlow/build/TwoPtFlow
export LD_LIBRARY_PATH=<lib_path>/OMNI/lib/:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
TotalCORES=${SLURM_NTASKS} 

echo "  ${MPIRUN} -n ${TotalCORES} ${MemOpt} ${HADRONS} $SLURM_SUBMIT_DIR/$input.xml --grid ${Lattice}  --mpi ${MPIGrid} --threads ${OMP_NUM_THREADS}  --comms-concurrent --comms-overlap  --shm 2048 --dslash-unroll --decomposition"
${MPIRUN} -n ${TotalCORES} ${MemOpt} ${HADRONS} $SLURM_SUBMIT_DIR/$input.xml --grid ${Lattice}  --mpi ${MPIGrid} --threads ${OMP_NUM_THREADS} --comms-concurrent --comms-overlap  --shm 2048 --dslash-unroll --decomposition

############################################################################################# 
############################################################################################# 


TotalToc=`date +%s`
echo "# time-finish "`date`

TotalTime=$(( $TotalToc - $TotalTic ))
TotalHours=`echo "$TotalTime / 3600" | bc -l`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Total time  $TotalTime [sec] = $TotalHours [h]"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
