#!/bin/bash --login
##
## SLURM job submission file example for running IFOS3D on a test simulation
## testcase. No profiling support available. 
##============================================================================

#SBATCH --ntasks=1000
#SBATCH --nodes=42
##SBATCH --partition=debugq
#SBATCH -t 02:30:00
#SBATCH -J ifos3d_xhole
#SBATCH --account=pawsey0001
#SBATCH --export=NONE

module load cray-hdf5-parallel

# run forward simulation to obtain observed data
make clean
make install 
srun -n 1000 -N 42 --export=ALL ../bin/ifos3d ./in_and_out/xhole_FW.json | tee ./in_and_out/xhole_FW.out

