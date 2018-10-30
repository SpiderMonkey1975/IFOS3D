#!/bin/bash --login
##
## SLURM job submission file example for running IFOS2D on a small simulation
## testcase. No profiling support available. 
##============================================================================

#SBATCH --ntasks=32
#SBATCH --nodes=2
#SBATCH --partition=debugq
#SBATCH -t 00:10:00
#SBATCH -J ifos3d_test
#SBATCH --account=pawsey0001
#SBATCH --export=NONE

#module load cray-hdf5-parallel perftools-base perftools
module load cray-hdf5-parallel 

# run forward simulation to obtain observed data
make clean
make install  MODEL=hh_toy_true.c
srun -n 32 -N 2 --export=ALL ../bin/ifos3d ./in_and_out/ifos3d_toy_FW.json | tee ./in_and_out/ifos3D.out
#srun -n 32 -N 2 --export=ALL ../bin/ifos3d+pat ./in_and_out/ifos3d_toy_FW.json | tee ./in_and_out/ifos3D.out

#cp model/toy.vs_it0 model/toy.vs.true
#cp model/toy.vp_it0 model/toy.vp.true
#cp model/toy.rho_it0 model/toy.rho.true


# invert observed data using homogeneous starting model
#make clean
#make install MODEL=hh_toy_start.c
#srun -n 32 -N 2 --export=ALL ../bin/ifos3d ./in_and_out/ifos3d_toy.json | tee ./in_and_out/ifos3D_INV.out
