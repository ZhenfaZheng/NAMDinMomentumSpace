#!/bin/bash
#SBATCH --partition          normal,normal2,normal3,normal4
#SBATCH --time               24:00:00
#SBATCH --nodes              1
#SBATCH --ntasks-per-node    28

# Load Modules
module load hdf5/1.12.1/intel
source /data/soft/intel/oneapi2022/setvars.sh
export PATH=/data/app/qe-7.2/Hefei-NAMD/NAMD-EPC/src:$PATH

echo "Job Running ..."

JOB='namd'
EXE='namd-epc'
mpirun -n $SLURM_NTASKS $EXE > $JOB.out

echo "Job Done!"
